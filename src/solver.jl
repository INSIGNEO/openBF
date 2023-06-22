"""
    pressure(A :: Float64, A0 :: Float64, beta :: Float64, Pext :: Float64)

Compute pressure at a specific node by means of the linear-elastic constitutive
equation.
"""
# function pressure(A :: Float64, A0 :: Float64, beta :: Float64, Pext :: Float64)
#     return Pext + beta*(sqrt(A/A0) - 1.0)
# end

pressure(A::Float64, A0::Float64, β::Float64, Pext::Float64) = Pext + β*(sqrt(A/A0) - 1.0) 


#muladd(β, √(A/A0)-1.0, Pext)

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
function muscl_(v :: Vessel, dt :: Float64, b :: Blood)
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

    v.uStarA[2:v.M+1] =view(v.vA,2:v.M+1).-invDxDt.*diff(v.fluxA)
    v.uStarQ[2:v.M+1] =view(v.vQ,2:v.M+1).-invDxDt.*diff(v.fluxQ)
    v.uStarA[1], v.uStarA[end] = v.uStarA[2], v.uStarA[end-1]
    v.uStarQ[1], v.uStarQ[end] = v.uStarQ[2], v.uStarQ[end-1]

    computeLimiter!(v.slopesA, v, v.uStarA, v.invDx, v.dU)
    computeLimiter!(v.slopesQ, v, v.uStarQ, v.invDx, v.dU)

    v.Al = v.uStarA.+v.slopesA.*v.halfDx
    v.Ar = v.uStarA.-v.slopesA.*v.halfDx
    v.Ql = v.uStarQ.+v.slopesQ.*v.halfDx
    v.Qr = v.uStarQ.-v.slopesQ.*v.halfDx

    computeFlux!(v.Fl, v.gamma_ghost, v.Al, v.Ql)
    computeFlux!(v.Fr, v.gamma_ghost, v.Ar, v.Qr)

    v.fluxA = (v.Qr[2:v.M+2].+v.Ql[1:v.M+1].-dxDt.*(v.Ar[2:v.M+2].-v.Al[1:v.M+1])).*0.5
    v.fluxQ = (v.Fr[2:v.M+2].+v.Fl[1:v.M+1].-dxDt.*(v.Qr[2:v.M+2].-v.Ql[1:v.M+1])).*0.5
    
    v.A = (v.A.+v.uStarA[2:v.M+1].-invDxDt.*diff(v.fluxA)).*0.5
    v.Q = (v.Q.+v.uStarQ[2:v.M+1].-invDxDt.*diff(v.fluxQ)).*0.5

    for i = 1:v.M

        v.P[i] = pressure(v.A[i], v.A0[i], v.beta[i], v.Pext)
        v.u[i] = v.Q[i]/v.A[i]
        v.c[i] = waveSpeed(v.A[i], v.gamma[i])
      end

    #source term
    dt_rho_inv = dt/b.rho

    for i = eachindex(v.A)
      s_A0_inv = sqrt(1.0/v.A0[i])
      s_A = sqrt(v.A[i])
      s_A_inv = 1.0/s_A
      A0_inv = 1/v.A0[i]
      # v.Q[i] += dt_rho_inv*( -v.viscT*v.Q[i]/v.A[i] +
      #           v.A[i]*(v.beta[i]*0.5*v.dA0dx[i] -
      #           (s_A*s_A0_inv-1.)*v.dTaudx[i]))
    # v.Q[i] -= v.viscT*v.Q[i]/(v.A[i])*dt_rho_inv   #viscosity
    #     # v.Q[i] += dt*0.5*v.beta[i]*v.A[i]^1.5/(v.A0[i]*b.rho)*v.dA0dx[i]   #dP/dA0
    # v.Q[i] += 0.5*v.beta[i]*v.A[i]^1.5*A0_inv*dt_rho_inv*v.dA0dx[i]   #dP/dA0
    # v.Q[i] -= (v.A[i]*dt_rho_inv)*(sqrt(v.A[i]*A0_inv)-1.)*v.dTaudx[i]  #dP/dh0
    #     # v.Q[i] += dt*v.A[i]*s_A*rho_inv *
    #     #           (0.5*v.beta[i]*s_A0_inv*s_A0_inv*v.dA0dx[i] -
    #     #           (s_A0_inv - s_A_inv)*v.dTaudx[i] -
    #     #            2*(b.gamma_profile + 2)*pi*b.mu*v.Q[i]*
    #     #            s_A_inv*s_A_inv*s_A_inv*s_A_inv*s_A_inv )

    # v.Q[i] += dt_rho_inv*( -v.viscT*v.Q[i]*s_A_inv*s_A_inv +
    #                     v.A[i]*(v.beta[i]*0.5*v.dA0dx[i] -
    #                             (s_A*s_A0_inv-1.)*v.dTaudx[i]))


        v.Q[i] -= v.viscT*v.Q[i]/(v.A[i]*b.rho)*dt   #viscosity
        v.Q[i] += dt*0.5*v.beta[i]*v.A[i]^1.5/(v.A0[i]*b.rho)*v.dA0dx[i]   #dP/dA0
        v.Q[i] -= dt*(v.A[i]/b.rho)*(sqrt(v.A[i]/v.A0[i])-1.)*v.dTaudx[i]  #dP/dh0

    end

    v.P = pressure.(v.A, v.A0, v.beta, v.Pext)
    v.c = waveSpeed.(v.A, v.gamma)
    

    #parabolic system (viscoelastic part)
    # if v.wallVa[1] != 0.0
    #     Td = 1.0/dt + v.wallVb
    #     Tlu = -v.wallVa
    #     T = Tridiagonal(Tlu[1:end-1], Td, Tlu[2:end])

    #     d = (1.0/dt - v.wallVb).*v.Q
    #     d[1] += v.wallVa[2]*v.Q[2]
    #     for i = 2:v.M-1
    #         d[i] += v.wallVa[i+1]*v.Q[i+1] + v.wallVa[i-1]*v.Q[i-1]
    #     end
    #     d[end] += v.wallVa[end-1]*v.Q[end-1]

    #     v.Q = T\d
    # end

    v.u = v.Q./v.A
end

function muscl(v :: Vessel, dt :: Float64, b :: Blood)

  v.vA[1] = v.U00A
  v.vA[end] = v.UM1A

  v.vQ[1] = v.U00Q
  v.vQ[end] = v.UM1Q

  for i = 2:v.M+1
    v.vA[i] = v.A[i-1]
    v.vQ[i] = v.Q[i-1]
  end

  v.slopesA = computeLimiter(v, v.vA, v.invDx, v.dU, v.slopesA)
  v.slopesQ = computeLimiter(v, v.vQ, v.invDx, v.dU, v.slopesQ)

  for i = 1:v.M+2
    v.Al[i] = v.vA[i] + v.slopesA[i]*v.halfDx
    v.Ar[i] = v.vA[i] - v.slopesA[i]*v.halfDx

    v.Ql[i] = v.vQ[i] + v.slopesQ[i]*v.halfDx
    v.Qr[i] = v.vQ[i] - v.slopesQ[i]*v.halfDx
  end

  v.Fl = computeFlux(v, v.Al, v.Ql, v.Fl)
  v.Fr = computeFlux(v, v.Ar, v.Qr, v.Fr)

  dxDt = v.dx/dt
  invDxDt = 1.0/dxDt

  for i = 1:v.M+1
    v.flux[1, i] = 0.5 * ( v.Fr[1, i+1] + v.Fl[1, i] - dxDt *
                        (v.Ar[i+1]    - v.Al[i]) )
    v.flux[2, i] = 0.5 * ( v.Fr[2, i+1] + v.Fl[2, i] - dxDt *
                        (v.Qr[i+1]    - v.Ql[i]) )
  end

  for i = 2:v.M+1
    v.uStar[1, i] = v.vA[i] + invDxDt * (v.flux[1, i-1] - v.flux[1, i])
    v.uStar[2, i] = v.vQ[i] + invDxDt * (v.flux[2, i-1] - v.flux[2, i])
  end

  v.uStar[1,1] = v.uStar[1,2]
  v.uStar[2,1] = v.uStar[2,2]
  v.uStar[1,end] = v.uStar[1,end-1]
  v.uStar[2,end] = v.uStar[2,end-1]

  v.slopesA = computeLimiter(v, vec(v.uStar[1,:]), v.invDx, v.dU, v.slopesA)
  v.slopesQ = computeLimiter(v, vec(v.uStar[2,:]), v.invDx, v.dU, v.slopesQ)

  for i = 1:v.M+2
    v.Al[i] = v.uStar[1,i] + v.slopesA[i]*v.halfDx
    v.Ar[i] = v.uStar[1,i] - v.slopesA[i]*v.halfDx

    v.Ql[i] = v.uStar[2,i] + v.slopesQ[i]*v.halfDx
    v.Qr[i] = v.uStar[2,i] - v.slopesQ[i]*v.halfDx
  end

  v.Fl = computeFlux(v, v.Al, v.Ql, v.Fl)
  v.Fr = computeFlux(v, v.Ar, v.Qr, v.Fr)

  for i = 1:v.M+1
    v.flux[1, i] = 0.5 * ( v.Fr[1, i+1] + v.Fl[1, i] - dxDt *
                        (v.Ar[i+1]    - v.Al[i]) )
    v.flux[2, i] = 0.5 * ( v.Fr[2, i+1] + v.Fl[2, i] - dxDt *
                        (v.Qr[i+1]    - v.Ql[i]) )
  end

  for i = 2:v.M+1
    v.A[i-1] = 0.5*(v.A[i-1] + v.uStar[1, i] + invDxDt *
                    (v.flux[1, i-1] - v.flux[1, i]) )
    v.Q[i-1] = 0.5*(v.Q[i-1] + v.uStar[2, i] + invDxDt *
                    (v.flux[2, i-1] - v.flux[2, i]))
  end

  #source term
  gamma_profile = 2
  for i = 1:v.M
    v.Q[i] -= 2*(gamma_profile+2)*pi*b.mu*v.Q[i]/(v.A[i]*b.rho)*dt   #viscosity
    v.Q[i] += dt*0.5*v.beta[i]*v.A[i]^1.5/(v.A0[i]*b.rho)*v.dA0dx[i]   #dP/dA0
    v.Q[i] -= dt*(v.A[i]/b.rho)*(sqrt(v.A[i]/v.A0[i])-1.)*v.dTaudx[i]  #dP/dh0

    v.P[i] = pressure(v.A[i], v.A0[i], v.beta[i], v.Pext)
    v.u[i] = v.Q[i]/v.A[i]
    v.c[i] = waveSpeed(v.A[i], v.gamma[i])
  end

end

function computeFlux(v :: Vessel, A :: Array{Float64, 1},
                     Q :: Array{Float64, 1}, Flux :: Array{Float64, 2})

  for i in 1:v.M+2
    Flux[1,i] = Q[i]
    Flux[2,i] = Q[i]*Q[i]/A[i] + v.gamma_ghost[i] * A[i]^1.5
  end

  return Flux
end

function computeLimiter(v :: Vessel, U :: Array{Float64, 1},
                        invDx :: Float64, dU :: Array{Float64, 2},
                        slopes :: Array{Float64, 1})

  for i = 2:v.M+2
    dU[1, i]   = (U[i] - U[i-1])*invDx
    dU[2, i-1] = (U[i] - U[i-1])*invDx
  end

  return superBee(v, dU, slopes)

end

function superBee(v :: Vessel, dU :: Array{Float64, 2}, slopes :: Array{Float64, 1})

  for i in 1:v.M+2
    s1 = minmod(dU[1,i], 2*dU[2,i])
    s2 = minmod(2*dU[1,i], dU[2,i])
    slopes[i] = maxmod(s1, s2)
  end

  return slopes
end

function maxmod(a :: Float64, b :: Float64)
  if abs(a) > abs(b)
    return a

  else
    return b

  end
end

function minmod(a :: Float64, b :: Float64)
  if a*b <= 0.
    return 0.

  elseif abs(a) < abs(b)
    return a

  else
    return b

  end
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
# maxMod(a::Float64, b::Float64) = max(a, b)

maxMod(a::Float64, b::Float64) = 0.5 * sum(sign, (a, b)) * maximum(abs.((a,b)))

"""
    minMod(a :: Float64, b :: Float64)
"""
# minMod(a::Float64, b::Float64) = a≤0.0 || b≤0.0 ? 0.0 : min(a, b)


minMod(a::Float64, b::Float64) = 0.5 * sum(sign, (a, b)) * minimum(abs.((a,b)))

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
        dU[1,i] = (U[i] - U[i-1])*v.invDx
    end
    for i = 1:v.M+1
        dU[2,i] = (U[i+1] - U[i])*v.invDx
    end
    slopes = superBee(dU)
end
