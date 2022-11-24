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
# function waveSpeed(A :: Float64, gamma :: Float64)
#     return sqrt(3*gamma*sqrt(A)*0.5)
# end

waveSpeed(A::Float64, γ::Float64) = √(1.5γ√A)

# function waveSpeedSA(sA :: Float64, gamma :: Float64)
#     return sqrt(1.5*gamma*sA)
# end


"""
    calculateDeltaT(vessels :: Array{Vessel,1}, Ccfl :: Float64)

Compute global minimum Δt through CFL condition. For each vessel the Δt is computed as

Delta t = C_{CFL} frac{Delta x}{S_{max}}

where Smax is the maximum between the forward and the backward characteristics, and Ccfl is the Courant-Friedrichs-Lewy condition defined by the user.
"""
calculateDeltaT(vessels::Array{Vessel,1}, C::Float64) = minimum([v.dx*C/maximum(abs.(v.u+v.c)) for v=vessels])


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

    v.slopesA = computeLimiter(v, v.vA, v.invDx, v.dU)
    v.slopesQ = computeLimiter(v, v.vQ, v.invDx, v.dU)

    v.Al = muladd(v.slopesA, v.halfDx, v.vA) 
    v.Ar = muladd(-v.slopesA, v.halfDx, v.vA) 
    v.Ql = muladd(v.slopesQ, v.halfDx, v.vQ) 
    v.Qr = muladd(-v.slopesQ, v.halfDx, v.vQ) 

    v.Fl = computeFlux(v, v.Al, v.Ql, v.Fl)
    v.Fr = computeFlux(v, v.Ar, v.Qr, v.Fr)

    dxDt = v.dx/dt
    invDxDt = 1.0/dxDt

    for i = 1:v.M+1
        v.flux[1,i] = 0.5*(v.Fr[1,i+1] + v.Fl[1,i] - dxDt*(v.Ar[i+1] - v.Al[i]))
        v.flux[2,i] = 0.5*(v.Fr[2,i+1] + v.Fl[2,i] - dxDt*(v.Qr[i+1] - v.Ql[i]))
    end

    for i = 2:v.M+1
        v.uStar[1,i] = v.vA[i] + invDxDt*(v.flux[1,i-1] - v.flux[1,i])
        v.uStar[2,i] = v.vQ[i] + invDxDt*(v.flux[2,i-1] - v.flux[2,i])
    end

    v.uStar[1,1] = v.uStar[1,2]
    v.uStar[2,1] = v.uStar[2,2]
    v.uStar[1,end] = v.uStar[1,end-1]
    v.uStar[2,end] = v.uStar[2,end-1]

    v.slopesA = computeLimiter(v, v.uStar, 1, v.invDx, v.dU)
    v.slopesQ = computeLimiter(v, v.uStar, 2, v.invDx, v.dU)

    for i = 1:v.M+2
        v.Al[i] = v.uStar[1,i] + v.slopesA[i]*v.halfDx
        v.Ar[i] = v.uStar[1,i] - v.slopesA[i]*v.halfDx

        v.Ql[i] = v.uStar[2,i] + v.slopesQ[i]*v.halfDx
        v.Qr[i] = v.uStar[2,i] - v.slopesQ[i]*v.halfDx
    end

    v.Fl = computeFlux(v, v.Al, v.Ql, v.Fl)
    v.Fr = computeFlux(v, v.Ar, v.Qr, v.Fr)

    for i = 1:v.M+1
        v.flux[1,i] = 0.5*(v.Fr[1,i+1] + v.Fl[1,i] - dxDt*(v.Ar[i+1] - v.Al[i]))
        v.flux[2,i] = 0.5*(v.Fr[2,i+1] + v.Fl[2,i] - dxDt*(v.Qr[i+1] - v.Ql[i]))
    end

    for i = 2:v.M+1
        v.A[i-1] = 0.5*(v.A[i-1] + v.uStar[1,i] + invDxDt*(v.flux[1,i-1] - v.flux[1,i]))
        v.Q[i-1] = 0.5*(v.Q[i-1] + v.uStar[2,i] + invDxDt*(v.flux[2,i-1] - v.flux[2,i]))
    end

    #source term
    for i = 1:v.M
        s_A = sqrt(v.A[i])
        Si = - v.viscT*v.Q[i]/v.A[i] - v.wallE[i]*(s_A - v.s_A0[i])*v.A[i]
        v.Q[i] += dt*Si

        v.P[i] = pressure(v.A[i], v.A0[1], v.beta[i], v.Pext)
        v.c[i] = waveSpeed(v.A[i], v.gamma[i])

    end

    #parabolic system (viscoelastic part)
    if v.wallVa[1] != 0.0
        Td = 1.0/dt + v.wallVb
        Tlu = -v.wallVa
        T = Tridiagonal(Tlu[1:end-1], Td, Tlu[2:end])

        d = (1.0/dt - v.wallVb).*v.Q
        d[1] += v.wallVa[2]*v.Q[2]
        for i = 2:v.M-1
            d[i] += v.wallVa[i+1]*v.Q[i+1] + v.wallVa[i-1]*v.Q[i-1]
        end
        d[end] += v.wallVa[end-1]*v.Q[end-1]

        v.Q = T\d
    end

    v.u = v.Q./v.A
end


"""
    computeFlux(v :: Vessel, A :: Array{Float64,1}, Q :: Array{Float64,1},
                Flux :: Array{Float64,2})
"""
function computeFlux(v :: Vessel, A :: Array{Float64,1}, Q :: Array{Float64,1},
                     Flux :: Array{Float64,2})
    Flux[1,:] = Q
    Flux[2,:] = Q.*Q./A .+ v.gamma_ghost.*A.*.√A
    Flux
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
    computeLimiter(v :: Vessel, U :: Array{Float64,1}, invDx :: Float64,
                   dU :: Array{Float64,2})
"""
function computeLimiter(v :: Vessel, U :: Array{Float64,1}, invDx :: Float64,
                        dU :: Array{Float64,2})
    for i = 2:v.M+2
        dU[1,i]   = (U[i] - U[i-1])*v.invDx
        dU[2,i-1] = dU[1,i]
    end
    superBee(dU)
end


function computeLimiter(v :: Vessel, U :: Array{Float64,2}, idx :: Int,
        invDx :: Float64, dU :: Array{Float64,2})
    U = U[idx,:]
    for i = 2:v.M+2
        dU[1,i]   = (U[i] - U[i-1])*invDx
        dU[2,i-1] = dU[1,i]
    end
    superBee(dU)
end
