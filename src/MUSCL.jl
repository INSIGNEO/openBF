# Documentation to be completed.
# See references[^1][^2][^3][^4]
function computeFlux(v :: Vessel, A :: Array{Float64, 1},
                     Q :: Array{Float64, 1}, Flux :: Array{Float64, 2})

  for i in 1:v.M+2
    Flux[1,i] = Q[i]
    Flux[2,i] = Q[i]*Q[i]/A[i] + v.gamma_ghost[i] * A[i]^1.5
  end

  return Flux
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

function superBee(v :: Vessel, dU :: Array{Float64, 2}, slopes :: Array{Float64, 1})

  for i in 1:v.M+2
    s1 = minmod(dU[1,i], 2*dU[2,i])
    s2 = minmod(2*dU[1,i], dU[2,i])
    slopes[i] = maxmod(s1, s2)
  end

  return slopes
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

function MUSCL(v :: Vessel, dt :: Float64, b :: Blood)

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
  invDxDt = 1./dxDt

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
  for i = 1:v.M
    v.Q[i] -= 2*(b.gamma_profile+2)*pi*b.mu*v.Q[i]/(v.A[i]*b.rho)*dt   #viscosity
    v.Q[i] += dt*0.5*v.beta[i]*v.A[i]^1.5/(v.A0[i]*b.rho)*v.dA0dx[i]   #dP/dA0
    v.Q[i] -= dt*(v.A[i]/b.rho)*(sqrt(v.A[i]/v.A0[i])-1.)*v.dTaudx[i]  #dP/dh0

    v.P[i] = pressure(v.A[i], v.A0[i], v.beta[i], v.Pext)
    v.u[i] = v.Q[i]/v.A[i]
    v.c[i] = waveSpeed(v.A[i], v.gamma[i])
  end

end
#
# ### References
#
# [^1]: LeVeque RJ. Finite volume methods for hyperbolic problems. Cambridge university press; 2002.
# [^2]: LeVeque RJ. Numerical methods for conservation laws. Birkhauser; 1992.
# [^3]: Toro, Eleuterio F. Riemann solvers and numerical methods for fluid dynamics: a practical introduction. Springer Science and Business Media, 2009.
# [^4]: Van Leer B. Towards the ultimate conservative difference scheme. Journal of Computational Physics. 1997;135(2):229-48.
