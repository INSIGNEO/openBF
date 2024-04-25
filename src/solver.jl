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


function calculateΔt(n::Network)
    minΔt = 1.0
    Δt = 1.0
    for (_, v) in n.vessels
        maxspeed = 0.0
        for i=eachindex(v.u)
            @inbounds speed = abs(v.u[i] + wave_speed(v.A[i], v.gamma[i+1]))
            maxspeed = ifelse(speed>maxspeed, speed, maxspeed)
        end
        Δt = v.dx / maxspeed
        minΔt = ifelse(Δt < minΔt, Δt, minΔt)
    end
    minΔt*n.Ccfl
end

function solve!(n::Network, dt::Float64, current_time::Float64)::Nothing
    for edge in n.edges
        s = Graphs.src(edge)
        t = Graphs.dst(edge)

        # inlet
        s == 1 && inbc!(n.vessels[(s, t)], current_time, dt, n.heart)

        # TODO: multiple inlets

        # vessel
        muscl!(n.vessels[(s, t)], dt, n.blood)

        # downstream
        outdeg::Int64 = Graphs.outdegree(n.graph, t)
        if outdeg == 0 # outlet
            outbc!(n.vessels[(s, t)], dt, n.blood.rho)
        elseif outdeg == 1
            indeg::Int64 = Graphs.indegree(n.graph, t)
            d::Int64 = first(Graphs.outneighbors(n.graph, t))
            if indeg == 1 # conjunction
                join_vessels!(n.vessels[(s, t)], n.vessels[(t, d)], n.blood.rho)
            
            # TODO: test ANASTOMOSIS!!!
            elseif indeg == 2 # anastomosis
                ps::Vector{Int64} = Graphs.inneighbors(n.graph, t)
                if t == max(ps[1], ps[2])
                    solveAnastomosis(
                        n.vessels[(ps[1], t)],
                        n.vessels[(ps[2], t)],
                        n.vessels[(t, d)],
                    )
                end
            end
            
        elseif outdeg == 2 # bifurcation
            ds::Vector{Int64} = Graphs.outneighbors(n.graph, t)
            join_vessels!(n.vessels[(s, t)], n.vessels[t, ds[1]], n.vessels[t, ds[2]])
        end
    end
end

function muscl!(v::Vessel, dt::Float64, b::Blood)
    dxDt = v.dx / dt
    invDxDt = 1.0 / dxDt

    # step-1
    # i = 1
    @inbounds v.vA[1] = v.U00A
    @inbounds v.vQ[1] = v.U00Q
    # 2:v.M+1
    for i = eachindex(v.A, v.Q)
        @inbounds v.vA[i+1] = v.A[i]
        @inbounds v.vQ[i+1] = v.Q[i]
    end
    # i = M+2
    @inbounds v.vA[end] = v.UM1A
    @inbounds v.vQ[end] = v.UM1Q

    #
    limiter!(v, v.vA, v.slopesA)
    limiter!(v, v.vQ, v.slopesQ)

    #
    for i = eachindex(v.Al)
        @inbounds v.Al[i] = v.vA[i] + v.slopesA[i]
        @inbounds v.Ar[i] = v.vA[i] - v.slopesA[i]

        @inbounds v.Ql[i] = v.vQ[i] + v.slopesQ[i]
        @inbounds v.Qr[i] = v.vQ[i] - v.slopesQ[i]

        @inbounds v.Fl[i] = v.Ql[i] * v.Ql[i] / v.Al[i] + v.gamma[i] * v.Al[i] * sqrt(v.Al[i])
        @inbounds v.Fr[i] = v.Qr[i] * v.Qr[i] / v.Ar[i] + v.gamma[i] * v.Ar[i] * sqrt(v.Ar[i])
    end

    #
    for i = 1:v.M+1
        @inbounds v.fluxA[i] = 0.5 * (v.Qr[i+1] + v.Ql[i] - dxDt * (v.Ar[i+1] - v.Al[i]))
        @inbounds v.fluxQ[i] = 0.5 * (v.Fr[i+1] + v.Fl[i] - dxDt * (v.Qr[i+1] - v.Ql[i]))
    end

    # step-2
    for i = 2:v.M+1
        @inbounds v.vA[i] += invDxDt * (v.fluxA[i-1] - v.fluxA[i])
        @inbounds v.vQ[i] += invDxDt * (v.fluxQ[i-1] - v.fluxQ[i])
    end
    # i=1
    @inbounds v.vA[1] = v.vA[2]
    @inbounds v.vQ[1] = v.vQ[2]
    # i=M+2
    @inbounds v.vA[end] = v.vA[end-1]
    @inbounds v.vQ[end] = v.vQ[end-1]

    #
    limiter!(v, v.vA, v.slopesA)
    limiter!(v, v.vQ, v.slopesQ)

    #
    for i = eachindex(v.Al)
        @inbounds v.Al[i] = v.vA[i] + v.slopesA[i]
        @inbounds v.Ar[i] = v.vA[i] - v.slopesA[i]

        @inbounds v.Ql[i] = v.vQ[i] + v.slopesQ[i]
        @inbounds v.Qr[i] = v.vQ[i] - v.slopesQ[i]

        @inbounds v.Fl[i] = v.Ql[i] * v.Ql[i] / v.Al[i] + v.gamma[i] * v.Al[i] * sqrt(v.Al[i])
        @inbounds v.Fr[i] = v.Qr[i] * v.Qr[i] / v.Ar[i] + v.gamma[i] * v.Ar[i] * sqrt(v.Ar[i])
    end

    #
    for i = 1:v.M+1
        @inbounds v.fluxA[i] = 0.5 * (v.Qr[i+1] + v.Ql[i] - dxDt * (v.Ar[i+1] - v.Al[i]))
        @inbounds v.fluxQ[i] = 0.5 * (v.Fr[i+1] + v.Fl[i] - dxDt * (v.Qr[i+1] - v.Ql[i]))
    end

    # 2:v.M+1
    for i = eachindex(v.A)
        @inbounds v.A[i] =
            0.5 * (v.A[i] + v.vA[i+1] + invDxDt * (v.fluxA[i] - v.fluxA[i+1]))
        @inbounds v.Q[i] =
            0.5 * (v.Q[i] + v.vQ[i+1] + invDxDt * (v.fluxQ[i] - v.fluxQ[i+1]))
    end

    #source term
    for i = eachindex(v.Q) # 1:v.M
        #viscosity
        @inbounds v.Q[i] -= 2 * (v.gamma_profile + 2) * pi * b.mu * v.Q[i] / (v.A[i] * b.rho) * dt

        if v.tapered
            #dP/dA0   
            @inbounds v.Q[i] += dt * 0.5 * v.beta[i] * sqrt(v.A[i])*v.A[i] / (v.A0[i] * b.rho) * v.dA0dx[i]

            #dP/dh0
            @inbounds v.Q[i] -= dt * (v.A[i] / b.rho) * (sqrt(v.A[i] / v.A0[i]) - 1.0) * v.dTaudx[i]
        end

        # update
        @inbounds v.u[i] = v.Q[i] / v.A[i]
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
end

function limiter!(
    v::Vessel,
    U::Vector{Float64},
    slopes::Vector{Float64},
)
    for i = 2:v.M+2
        @inbounds v.dUA[i] = (U[i] - U[i-1]) * v.invDx
        @inbounds v.dUQ[i-1] = (U[i] - U[i-1]) * v.invDx
    end
    superbee!(slopes, v.dUA, v.dUQ, v.halfDx)
end

# https://discourse.julialang.org/t/optimising-superbee-function/112568/12
function superbee!(s::Vector{Float64}, a::Vector{Float64}, b::Vector{Float64}, h::Float64)
    @simd for i = eachindex(s)
        # @inbounds is inferred automatically - yay for safety AND speed!
        ai = a[i]
        bi = b[i]
        t1 = max(min(ai,2bi),min(2ai,bi))
        t2 = min(max(ai,2bi),max(2ai,bi))
        s[i] = ifelse(ai>0, ifelse(bi>0, t1, zero(Float64)), ifelse(bi<0, t2, zero(Float64)))*h
    end
end
