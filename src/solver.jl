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
            set_outlet_bc(dt, n.vessels[(s, t)])
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
        v.vA[i] = v.A[i-1]
        v.vQ[i] = v.Q[i-1]
    end

    limiter!(v, v.vA, v.invDx, v.dU, v.slopesA)
    limiter!(v, v.vQ, v.invDx, v.dU, v.slopesQ)

    for i = 1:v.M+2
        v.Al[i] = v.vA[i] + v.slopesA[i] * v.halfDx
        v.Ar[i] = v.vA[i] - v.slopesA[i] * v.halfDx

        v.Ql[i] = v.vQ[i] + v.slopesQ[i] * v.halfDx
        v.Qr[i] = v.vQ[i] - v.slopesQ[i] * v.halfDx
    end

    flux!(v, v.Al, v.Ql, v.Fl)
    flux!(v, v.Ar, v.Qr, v.Fr)

    dxDt = v.dx / dt
    invDxDt = 1.0 / dxDt

    for i = 1:v.M+1
        v.flux[1, i] = 0.5 * (v.Fr[1, i+1] + v.Fl[1, i] - dxDt * (v.Ar[i+1] - v.Al[i]))
        v.flux[2, i] = 0.5 * (v.Fr[2, i+1] + v.Fl[2, i] - dxDt * (v.Qr[i+1] - v.Ql[i]))
    end

    for i = 2:v.M+1
        v.uStar[1, i] = v.vA[i] + invDxDt * (v.flux[1, i-1] - v.flux[1, i])
        v.uStar[2, i] = v.vQ[i] + invDxDt * (v.flux[2, i-1] - v.flux[2, i])
    end

    v.uStar[1, 1] = v.uStar[1, 2]
    v.uStar[2, 1] = v.uStar[2, 2]
    v.uStar[1, end] = v.uStar[1, end-1]
    v.uStar[2, end] = v.uStar[2, end-1]

    limiter!(v, vec(v.uStar[1, :]), v.invDx, v.dU, v.slopesA)
    limiter!(v, vec(v.uStar[2, :]), v.invDx, v.dU, v.slopesQ)

    for i = 1:v.M+2
        v.Al[i] = v.uStar[1, i] + v.slopesA[i] * v.halfDx
        v.Ar[i] = v.uStar[1, i] - v.slopesA[i] * v.halfDx

        v.Ql[i] = v.uStar[2, i] + v.slopesQ[i] * v.halfDx
        v.Qr[i] = v.uStar[2, i] - v.slopesQ[i] * v.halfDx
    end

    flux!(v, v.Al, v.Ql, v.Fl)
    flux!(v, v.Ar, v.Qr, v.Fr)

    for i = 1:v.M+1
        v.flux[1, i] = 0.5 * (v.Fr[1, i+1] + v.Fl[1, i] - dxDt * (v.Ar[i+1] - v.Al[i]))
        v.flux[2, i] = 0.5 * (v.Fr[2, i+1] + v.Fl[2, i] - dxDt * (v.Qr[i+1] - v.Ql[i]))
    end

    for i = 2:v.M+1
        v.A[i-1] =
            0.5 * (v.A[i-1] + v.uStar[1, i] + invDxDt * (v.flux[1, i-1] - v.flux[1, i]))
        v.Q[i-1] =
            0.5 * (v.Q[i-1] + v.uStar[2, i] + invDxDt * (v.flux[2, i-1] - v.flux[2, i]))
    end

    #source term
    for i = 1:v.M
        v.Q[i] -= 2 * (v.gamma_profile + 2) * pi * b.mu * v.Q[i] / (v.A[i] * b.rho) * dt   #viscosity
        v.Q[i] += dt * 0.5 * v.beta[i] * sqrt(v.A[i])*v.A[i] / (v.A0[i] * b.rho) * v.dA0dx[i]   #dP/dA0
        v.Q[i] -= dt * (v.A[i] / b.rho) * (sqrt(v.A[i] / v.A0[i]) - 1.0) * v.dTaudx[i]  #dP/dh0
    end

    #parabolic system (visco-elastic)
    if v.Cv[1] != 0.0
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

    for i=1:v.M
        v.P[i] = pressure(v.A[i], v.A0[i], v.beta[i], v.Pext)
        if (v.Cv[i] == 0) && (v.M != i != 1)
            v.P[i] -= v.Cv[i] * b.rho / v.A[i] * (v.Q[i] - v.Q[i-1])/v.dx
        end
        v.u[i] = v.Q[i] / v.A[i]
        v.c[i] = wave_speed(v.A[i], v.gamma[i])
    end
end

function limiter!(
    v::Vessel,
    U::Vector{Float64},
    invDx::Float64,
    dU::Array{Float64,2},
    slopes::Vector{Float64},
)
    for i = 2:v.M+2
        dU[1, i] = (U[i] - U[i-1]) * invDx
        dU[2, i-1] = (U[i] - U[i-1]) * invDx
    end
    superbee!(v, dU, slopes)
end

function flux!(v::Vessel, A::Vector{Float64}, Q::Vector{Float64}, Flux::Array{Float64,2})
    for i = 1:v.M+2
        Flux[1, i] = Q[i]
        Flux[2, i] = Q[i] * Q[i] / A[i] + v.gamma_ghost[i] * A[i]*sqrt(A[i])
    end
end

maxmod(v) = 0.5 * sum(sign, v) * maximum(abs, v)
minmod(v) = 0.5 * sum(sign, v) * minimum(abs, v)
function superbee!(v::Vessel, dU::Array{Float64,2}, slopes::Vector{Float64})
    for i = 1:v.M+2
        slopes[i] =
            maxmod((minmod((dU[1, i], 2 * dU[2, i])), minmod((2 * dU[1, i], dU[2, i]))))
    end
end
