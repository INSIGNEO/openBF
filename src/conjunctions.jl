#=
Copyright 2018 INSIGNEO Institute for in silico Medicine

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


function solveConjunction(b :: Blood, v1 :: Vessel, v2 :: Vessel)

  U = [v1.u[end],
       v2.u[ 1 ],
       sqrt(sqrt(v1.A[end])),
       sqrt(sqrt(v2.A[ 1 ]))]

  k1 = sqrt(0.5*3*v1.gamma[end])
  k2 = sqrt(0.5*3*v2.gamma[ 1 ])
  k  = [k1, k2]

  W = calculateWstarConj(U, k)
  J = calculateJacobianConj(b, v1, v2, U, k)
  F = calculateFofUconj(b, v1, v2, U, k, W)

  nr_toll_U = 1.e-5
  nr_toll_F = 1.e-5

  while true
    dU = J\(-F)
    U_new = U + dU

    if any(isnan(dot(F,F)))
      println(F)
      break
    end

    u_ok = 0
    f_ok = 0
    for i in 1:length(dU)
      if abs(dU[i]) <= nr_toll_U || abs(F[i]) <= nr_toll_F
        u_ok += 1
        f_ok += 1
      end
    end

    if u_ok == length(dU) || f_ok == length(dU)
      U = U_new
      break

    else
      U = U_new
      W = calculateWstarConj(U, k)
      F = calculateFofUconj(b, v1, v2, U, k, W)
    end
  end

  v1.u[end] = U[1]
  v2.u[ 1 ] = U[2]

  v1.A[end] = U[3]*U[3]*U[3]*U[3]
  v1.Q[end] = v1.u[end]*v1.A[end]

  v2.A[ 1 ] = U[4]^4
  v2.Q[ 1 ] = v2.u[1]*v2.A[1]

  v1.P[end] = pressure(v1.A[end], v1.A0[end], v1.beta[end], v1.Pext)
  v2.P[ 1 ] = pressure(v2.A[ 1 ], v2.A0[ 1 ], v2.beta[ 1 ], v2.Pext)

  v1.c[end] = waveSpeed(v1.A[end], v1.gamma[end])
  v2.c[ 1 ] = waveSpeed(v2.A[ 1 ], v2.gamma[ 1 ])
end


function calculateWstarConj(U :: Array, k :: Array)

  W1 = U[1] + 4*k[1] * U[3]
  W2 = U[2] - 4*k[2] * U[4]

  return [W1, W2]
end


function calculateFofUconj(b :: Blood, v1 :: Vessel, v2 :: Vessel,
                           U :: Array, k :: Array, W :: Array)

  f1 = U[1] + 4*k[1]*U[3] - W[1]

  f2 = U[2] - 4*k[2]*U[4] - W[2]

  f3 = U[1]*(U[3]*U[3]*U[3]*U[3]) - U[2]*(U[4]*U[4]*U[4]*U[4])

  f4 = 0.5*b.rho*U[1]*U[1] + v1.beta[end]*(U[3]*U[3]/sqrt(v1.A0[end]) - 1) -
     ( 0.5*b.rho*U[2]*U[2] + v2.beta[ 1 ]*(U[4]*U[4]/sqrt(v2.A0[ 1 ]) - 1) )

  return [f1, f2, f3, f4]
end

function calculateJacobianConj(b :: Blood, v1 :: Vessel, v2 :: Vessel,
                               U :: Array, k :: Array)

  J = eye(4)

  J[1,3] =  4*k[1]
  J[2,4] = -4*k[2]

  J[3,1] =  U[3]*U[3]*U[3]*U[3]
  J[3,2] = -U[4]*U[4]*U[4]*U[4]
  J[3,3] =  4*U[1]*(U[3]*U[3]*U[3]^3)
  J[3,4] = -4*U[2]*(U[4]*U[4]*U[4]^3)

  J[4,3] =  2*v1.beta[end]*U[3]/sqrt(v1.A0[end])
  J[4,4] = -2*v2.beta[ 1 ]*U[4]/sqrt(v2.A0[ 1 ])

  J[4,1] =  b.rho*U[1]
  J[4,2] = -b.rho*U[2]

  return J
end

# ### References
#
# [^1]: Press, William H. Numerical recipes 3rd edition:
# The art of scientific computing. Cambridge university press, 2007.
