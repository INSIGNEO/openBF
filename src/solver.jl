#=
Copyright 2022 INSIGNEO Institute for in silico Medicine

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


"""
    pressure(A :: Float64, A0 :: Float64, beta :: Float64, Pext :: Float64)

Compute pressure at a specific node by means of the linear-elastic constitutive
equation.
"""
# function pressure(A :: Float64, A0 :: Float64, beta :: Float64, Pext :: Float64)
#     return Pext + beta*(sqrt(A/A0) - 1.0)
# end

pressure(A::Float64, A0::Float64, β::Float64, Pext::Float64) = muladd(β, √(A/A0)-1.0, Pext)

# # version with pre-computed `sqrt(A/A0)`
# function pressure(s_A_over_A0 :: Float64, beta :: Float64, Pext :: Float64)
#     return Pext + beta*(s_A_over_A0 - 1.0)
# end


"""
    waveSpeed(A :: Float64, gamma :: Float64)

Compute pulse wave velocity at the given node.
"""
waveSpeed(A::Float64, γ::Float64) = √(1.5γ√A)

"""
    calculateDeltaT(vessels :: Array{Vessel,1}, Ccfl :: Float64)

Compute global minimum Δt through CFL condition. For each vessel the Δt is computed as

Delta t = C_{CFL} frac{Delta x}{S_{max}}

where Smax is the maximum between the forward and the backward characteristics, and Ccfl is the Courant-Friedrichs-Lewy condition defined by the user.
"""
calculateDeltaT(v::Vessel, C::Float64) = v.dx*C/maximum(abs.(v.u+v.c)) 
calculateDeltaT(v::Array{Vessel,1}, C::Float64) = minimum(calculateDeltaT.(v,C))



"""
    solveModel(vessel :: Array){Vessel,1}, edges :: Array{Int,2}, blood :: Blood,
               dt :: Float64, current_time :: Float64)

Run the solver on each vessel one-by-one as listed in the .yml file.
"""
function solveModel(vessels :: Array{Vessel,1}, edges :: Array{Int,2}, blood :: Blood,
                    dt :: Float64, current_time :: Float64)
    for j in 1:size(edges)[1]
        i = edges[j,1]
        v = vessels[i]

        solveVessel(v, blood, dt, current_time)
        solveOutlet(j, v, blood, vessels, edges, dt)
    end
end


"""
    solveVessel(vessel :: Vessel, blood :: Blood, dt :: Float64, current_time :: Float64)

Apply inlet boundary condition (if any) and run MUSCL solver along the vessel.
"""
function solveVessel(vessel :: Vessel, blood :: Blood, dt :: Float64,
                     current_time :: Float64)
    vessel.inlet && setInletBC(current_time, dt, vessel)
    muscl(vessel, dt, blood)
end


"""
    solveOutlet(j :: Int, vessel :: Vessel, blood :: Blood, vessels :: Array{Vessel,1},
                edges :: Array{Int,2})

Join vessel in case of junction or apply outlet boundary condition otherwise.
"""
function solveOutlet(j :: Int, vessel :: Vessel, blood :: Blood, vessels :: Array{Vessel,1},
                     edges :: Array{Int,2}, dt :: Float64)
    i = edges[j,1]
    t = edges[j,3]

    if vessel.outlet != "none"
        setOutletBC(dt, vessel)

    elseif size(findall(edges[:,2] .== t))[1] == 2
        d1_i = findall(edges[:,2] .== t)[1]
        d2_i = findall(edges[:,2] .== t)[2]
        joinVessels(blood, vessel, vessels[d1_i], vessels[d2_i])

    elseif size(findall(edges[:,3] .== t))[1] == 1
        d_i = findall(edges[:,2] .== t)[1]
        joinVessels(blood, vessel, vessels[d_i])

    elseif size(findall(edges[:,3] .== t))[1] == 2
        p1_i = findall(edges[:,3] .== t)[1]
        p2_i = findall(edges[:,3] .== t)[2]

        if maximum([p1_i, p2_i]) == i
            p2_i = minimum([p1_i, p2_i])
            d = findall(edges[:,2] .== t)[1]
            solveAnastomosis(vessel, vessels[p2_i], vessels[d])
        end
    end
end


"""
    muscl(v :: Vessel, dt :: Float64, b :: Blood)

Run MUSCL solver along the vessel.
"""
function muscl(v :: Vessel, dt :: Float64, b :: Blood)
    v.vA[1], v.vA[end] = v.U00A, v.UM1A
    v.vQ[1], v.vQ[end] = v.U00Q, v.UM1Q
    v.vA[2:v.M+1] = v.A
    v.vQ[2:v.M+1] = v.Q

    computeLimiter!(v.slopesA, v, v.vA, v.invDx, v.dU)
    computeLimiter!(v.slopesQ, v, v.vQ, v.invDx, v.dU)

    v.Al = muladd(v.slopesA, v.halfDx, v.vA) 
    v.Ar = muladd(-v.slopesA, v.halfDx, v.vA) 
    v.Ql = muladd(v.slopesQ, v.halfDx, v.vQ) 
    v.Qr = muladd(-v.slopesQ, v.halfDx, v.vQ) 

    computeFlux!(v.Fl, v.gamma_ghost, v.Al, v.Ql)
    computeFlux!(v.Fr, v.gamma_ghost, v.Ar, v.Qr)

    dxDt = v.dx/dt
    invDxDt = 1.0/dxDt

    v.fluxA = (v.Qr[2:v.M+2].+v.Ql[1:v.M+1].-dxDt.*(v.Ar[2:v.M+2].-v.Al[1:v.M+1])).*0.5
    v.fluxQ = (v.Fr[2:v.M+2].+v.Fl[1:v.M+1].-dxDt.*(v.Qr[2:v.M+2].-v.Ql[1:v.M+1])).*0.5

    # TODO: replace uStar matrix with two vectors
    v.uStar[1,2:v.M+1] = v.vA[2:v.M+1].-invDxDt.*diff(v.fluxA)
    v.uStar[2,2:v.M+1] = v.vQ[2:v.M+1].-invDxDt.*diff(v.fluxQ)
    v.uStar[1,1] = v.uStar[1,2]
    v.uStar[2,1] = v.uStar[2,2]
    v.uStar[1,end] = v.uStar[1,end-1]
    v.uStar[2,end] = v.uStar[2,end-1]

    computeLimiter!(v.slopesA, v, v.uStar[1,:], v.invDx, v.dU)
    computeLimiter!(v.slopesQ, v, v.uStar[2,:], v.invDx, v.dU)

    v.Al = v.uStar[1,:].+v.slopesA.*v.halfDx
    v.Ar = v.uStar[1,:].-v.slopesA.*v.halfDx
    v.Ql = v.uStar[2,:].+v.slopesQ.*v.halfDx
    v.Qr = v.uStar[2,:].-v.slopesQ.*v.halfDx

    computeFlux!(v.Fl, v.gamma_ghost, v.Al, v.Ql)
    computeFlux!(v.Fr, v.gamma_ghost, v.Ar, v.Qr)

    v.fluxA = (v.Qr[2:v.M+2].+v.Ql[1:v.M+1].-dxDt.*(v.Ar[2:v.M+2].-v.Al[1:v.M+1])).*0.5
    v.fluxQ = (v.Fr[2:v.M+2].+v.Fl[1:v.M+1].-dxDt.*(v.Qr[2:v.M+2].-v.Ql[1:v.M+1])).*0.5
    
    v.A = (v.A.+v.uStar[1,2:v.M+1].-invDxDt.*diff(v.fluxA)).*0.5
    v.Q = (v.Q.+v.uStar[2,2:v.M+1].-invDxDt.*diff(v.fluxQ)).*0.5

    #source term
    dt_rho_inv = dt/b.rho
    for i = eachindex(v.A)
      s_A0_inv = sqrt(1.0/v.A0[i])
      s_A = sqrt(v.A[i])
      s_A_inv = 1.0/s_A
      v.Q[i] += dt_rho_inv*( -v.viscT*v.Q[i]/v.A[i] +
                v.A[i]*(v.beta[i]*0.5*v.dA0dx[i] -
                (s_A*s_A0_inv-1.)*v.dTaudx[i]))
    end

    v.P = pressure.(v.A, v.A0, v.beta, v.Pext)
    v.c = waveSpeed.(v.A, v.gamma)
    

    ##parabolic system (viscoelastic part)
    #if v.wallVa[1] != 0.0
    #    Td = 1.0/dt + v.wallVb
    #    Tlu = -v.wallVa
    #    T = Tridiagonal(Tlu[1:end-1], Td, Tlu[2:end])

    #    d = (1.0/dt - v.wallVb).*v.Q
    #    d[1] += v.wallVa[2]*v.Q[2]
    #    for i = 2:v.M-1
    #        d[i] += v.wallVa[i+1]*v.Q[i+1] + v.wallVa[i-1]*v.Q[i-1]
    #    end
    #    d[end] += v.wallVa[end-1]*v.Q[end-1]

    #    v.Q = T\d
    #end

    v.u = v.Q./v.A
end


"""
    computeFlux(F::Vector{Float64}, gamma_ghost::Vector{Float64}, A::Vector{Float64}, Q::Vector{Float64})
"""
function computeFlux!(F::Vector{Float64}, gamma_ghost::Vector{Float64}, A::Vector{Float64}, Q::Vector{Float64})
    F[:] = Q.*Q./A.+gamma_ghost.*A.*.√A
end

"""
    maxMod(a :: Float64, b :: Float64)
"""
maxMod(a::Float64, b::Float64) = max(a, b)


"""
    minMod(a :: Float64, b :: Float64)
"""
minMod(a::Float64, b::Float64) = a≤0.0 || b≤0.0 ? 0.0 : min(a, b)


"""
    superBee(dU :: Array{Float64,2})
"""
function superBee(dU::Array{Float64,2})
    s1 = minMod.(dU[1,:], 2dU[2,:])
    s2 = minMod.(2dU[1,:], dU[2,:])
    maxMod.(s1, s2)
end


"""
    computeLimiter(v :: Vessel, U :: Vector{Float64}, invDx :: Float64,
                   dU :: Array{Float64,2})
"""
function computeLimiter!(slopes::Vector{Float64}, v :: Vessel, U :: Vector{Float64}, invDx :: Float64,
                        dU :: Array{Float64,2})
    for i = 2:v.M+2
        dU[1,i]   = (U[i] - U[i-1])*v.invDx
        dU[2,i-1] = dU[1,i]
    end
    slopes = superBee(dU)
end
