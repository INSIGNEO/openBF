#=
Copyright 2015-2024 INSIGNEO Institute for in silico Medicine

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=#

calculate_Δt(n::Network) =
    minimum(v.dx * v.Ccfl / maximum(abs, v.u .+ v.c) for (_, v) in n.vessels)

function solve!(n::Network, dt::Float64, current_time::Float64)
    for edge in edges(n.graph)
        s = src(edge)
        t = dst(edge)

        # inlet
        s == 1 && set_inlet_bc(current_time, dt, n.vessels[(s, t)], n.heart)

        # TODO: multiple inlets

        # vessel
        muscl!(n.vessels[(s, t)], dt, n.blood)

        if outdegree(n.graph, t) == 0 # outlet
            set_outlet_bc(dt, n.vessels[(s, t)], n.blood.rho)
        elseif outdegree(n.graph, t) == 1
            if indegree(n.graph, t) == 1 # conjunction
                d = first(outneighbors(n.graph, t))
                join_vessels!(n.blood, n.vessels[(s, t)], n.vessels[(t, d)])
                # TODO: test ANASTOMOSIS!!!
            elseif indegree(n.graph, t) == 2 # anastomosis
                p1, p2 = inneighbors(n.graph, t)
                if t == max(p1, p2)
                    d = outneighbors(grado, t)
                    solveAnastomosis(
                        n.vessels[(p1, t)],
                        n.vessels[(p2, t)],
                        n.vessels[(t, d)],
                    )
                end
            end
        elseif outdegree(n.graph, t) == 2 # bifurcation
            d1, d2 = outneighbors(n.graph, t)
            join_vessels!(n.vessels[(s, t)], n.vessels[t, d1], n.vessels[t, d2])
        end
    end
end

function muscl!(v::Vessel, dt::Float64, b::Blood)
    v.vA[1] = v.U00A
    v.vA[end] = v.UM1A

    v.vQ[1] = v.U00Q
    v.vQ[end] = v.UM1Q

    for i = 2:v.M+1
        @inbounds v.vA[i] = v.A[i-1]
        @inbounds v.vQ[i] = v.Q[i-1]
    end

    limiter!(v, v.vA, v.invDx, v.dUA, v.dUQ, v.slopesA)
    limiter!(v, v.vQ, v.invDx, v.dUA, v.dUQ, v.slopesQ)

    for i = eachindex(v.Al) # 1:v.M+2
        @inbounds v.Al[i] = v.vA[i] + v.slopesA[i]
        @inbounds v.Ar[i] = v.vA[i] - v.slopesA[i]

        @inbounds v.Ql[i] = v.vQ[i] + v.slopesQ[i]
        @inbounds v.Qr[i] = v.vQ[i] - v.slopesQ[i]

        @inbounds v.Fl[i] = v.Ql[i] * v.Ql[i] / v.Al[i] + v.gamma_ghost[i] * v.Al[i] * sqrt(v.Al[i])
        @inbounds v.Fr[i] = v.Qr[i] * v.Qr[i] / v.Ar[i] + v.gamma_ghost[i] * v.Ar[i] * sqrt(v.Ar[i])
    end

    dxDt = v.dx / dt
    invDxDt = 1.0 / dxDt

    for i = 1:v.M+1
        @inbounds v.fluxA[i] = 0.5 * (v.Qr[i+1] + v.Ql[i] - dxDt * (v.Ar[i+1] - v.Al[i]))
        @inbounds v.fluxQ[i] = 0.5 * (v.Fr[i+1] + v.Fl[i] - dxDt * (v.Qr[i+1] - v.Ql[i]))
    end

    for i = 2:v.M+1
        @inbounds v.uStarA[i] = v.vA[i] + invDxDt * (v.fluxA[i-1] - v.fluxA[i])
        @inbounds v.uStarQ[i] = v.vQ[i] + invDxDt * (v.fluxQ[i-1] - v.fluxQ[i])
    end

    v.uStarA[1] = v.uStarA[2]
    v.uStarQ[1] = v.uStarQ[2]
    v.uStarA[end] = v.uStarA[end-1]
    v.uStarQ[end] = v.uStarQ[end-1]

    limiter!(v, v.uStarA, v.invDx, v.dUA, v.dUQ, v.slopesA)
    limiter!(v, v.uStarQ, v.invDx, v.dUA, v.dUQ, v.slopesQ)

    for i = eachindex(v.Al) # 1:v.M+2
        @inbounds v.Al[i] = v.uStarA[i] + v.slopesA[i]
        @inbounds v.Ar[i] = v.uStarA[i] - v.slopesA[i]

        @inbounds v.Ql[i] = v.uStarQ[i] + v.slopesQ[i]
        @inbounds v.Qr[i] = v.uStarQ[i] - v.slopesQ[i]

        @inbounds v.Fl[i] = v.Ql[i] * v.Ql[i] / v.Al[i] + v.gamma_ghost[i] * v.Al[i] * sqrt(v.Al[i])
        @inbounds v.Fr[i] = v.Qr[i] * v.Qr[i] / v.Ar[i] + v.gamma_ghost[i] * v.Ar[i] * sqrt(v.Ar[i])
    end

    for i = 1:v.M+1
        @inbounds v.fluxA[i] = 0.5 * (v.Qr[i+1] + v.Ql[i] - dxDt * (v.Ar[i+1] - v.Al[i]))
        @inbounds v.fluxQ[i] = 0.5 * (v.Fr[i+1] + v.Fl[i] - dxDt * (v.Qr[i+1] - v.Ql[i]))
    end

    for i = 2:v.M+1
        @inbounds v.A[i-1] =
            0.5 * (v.A[i-1] + v.uStarA[i] + invDxDt * (v.fluxA[i-1] - v.fluxA[i]))
        @inbounds v.Q[i-1] =
            0.5 * (v.Q[i-1] + v.uStarQ[i] + invDxDt * (v.fluxQ[i-1] - v.fluxQ[i]))
    end

    #source term
    for i = eachindex(v.Q) # 1:v.M
        @inbounds v.Q[i] -= 2 * (v.gamma_profile + 2) * pi * b.mu * v.Q[i] / (v.A[i] * b.rho) * dt        #viscosity
        @inbounds v.Q[i] += dt * 0.5 * v.beta[i] * sqrt(v.A[i])*v.A[i] / (v.A0[i] * b.rho) * v.dA0dx[i]   #dP/dA0
        @inbounds v.Q[i] -= dt * (v.A[i] / b.rho) * (sqrt(v.A[i] / v.A0[i]) - 1.0) * v.dTaudx[i]          #dP/dh0
    end

    #parabolic system (visco-elastic)
    if v.viscoelastic
        a = v.Cv.*dt/(v.dx*v.dx)
        Tupper = -a[2:end]
        Tlower = -a[1:end-1]
        Tdiagonal = 1.0.+2.0.*a
        Tdiagonal[1] -= a[1]
        Tdiagonal[end] -= a[end]
        
        T = Tridiagonal(Tupper, Tdiagonal, Tlower)

        d = (1.0 .- 2.0.*a).*v.Q
        d[1] += a[2]*v.Q[2] + a[1]*v.Q[1]
        d[2:end-1] .+= a[1:end-2].*v.Q[1:end-2] .+ a[3:end].*v.Q[3:end]
        d[end] += a[end-1]*v.Q[end-1] + a[end]*v.Q[end]

        v.Q = T\d
    end

    #update
    for i=eachindex(v.u) # 1:v.M
        @inbounds v.u[i] = v.Q[i] / v.A[i]
        @inbounds v.c[i] = wave_speed(v.A[i], v.gamma[i])
    end
end

function limiter!(
    v::Vessel,
    U::Vector{Float64},
    invDx::Float64,
    dUA::Vector{Float64},
    dUQ::Vector{Float64},
    slopes::Vector{Float64},
)
    for i = 2:v.M+2
        @inbounds dUA[i] = (U[i] - U[i-1]) * invDx
        @inbounds dUQ[i-1] = (U[i] - U[i-1]) * invDx
    end
    superbee!(slopes, dUA, dUQ, v.halfDx)
end

maxmod(a::Float64, b::Float64) = 0.5*(sign(a) + sign(b))*max(abs(a), abs(b))
minmod(a::Float64, b::Float64) = 0.5*(sign(a) + sign(b))*min(abs(a), abs(b))
function superbee!(slopes::Vector{Float64}, dUA::Vector{Float64}, dUQ::Vector{Float64}, halfDx::Float64)
    for i = eachindex(slopes)
        @inbounds slopes[i] = maxmod(minmod(dUA[i], 2.0 * dUQ[i]), minmod(2.0 * dUA[i], dUQ[i])) * halfDx
    end
end
