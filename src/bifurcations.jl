# #=
# Copyright 2017 INSIGNEO Institute for in silico Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# =#
#
# function solveBifurcation(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel)
#
#   U = [v1.u[end],
#        v2.u[1],
#        v3.u[1],
#        sqrt(sqrt(v1.A[end])),
#        sqrt(sqrt(v2.A[1])),
#        sqrt(sqrt(v3.A[1]))]
#
#   k1 = sqrt(0.5*3*v1.gamma[end])
#   k2 = sqrt(0.5*3*v2.gamma[1])
#   k3 = sqrt(0.5*3*v3.gamma[1])
#   k = [k1, k2, k3]
#
#   W = calculateWstarBif(U, k)
#   J = calculateJacobianBif(v1, v2, v3, U, k)
#   F = calculateFofUBif(v1, v2, v3, U, k, W)
#
#   nr_toll_U = 1e-5
#   nr_toll_F = 1e-5
#
#   while true
#     dU = J\(-F)
#     U_new = U + dU
#
#     if any(isnan(dot(F,F)))
#       println(F)
#       @printf "error at bifurcation with vessels %s, %s, and %s \n" v1.label v2.label v3.label
#       break
#     end
#
#     u_ok = 0
#     f_ok = 0
#     for i in 1:length(dU)
#       if abs(dU[i]) <= nr_toll_U || abs(F[i]) <= nr_toll_F
#         u_ok += 1
#         f_ok += 1
#       end
#     end
#
#     if u_ok == length(dU) || f_ok == length(dU)
#       U = U_new
#       break
#     else
#       U = U_new
#       W = calculateWstarBif(U, k)
#       F = calculateFofUBif(v1, v2, v3, U, k, W)
#     end
#   end
#
#   v1.u[end] = U[1]
#   v2.u[1] = U[2]
#   v3.u[1] = U[3]
#
#   v1.A[end] = U[4]*U[4]*U[4]*U[4]
#   v2.A[1] = U[5]*U[5]*U[5]*U[5]
#   v3.A[1] = U[6]*U[6]*U[6]*U[6]
#
#   v1.Q[end] = v1.u[end]*v1.A[end]
#   v2.Q[1] = v2.u[1]*v2.A[1]
#   v3.Q[1] = v3.u[1]*v3.A[1]
#
#   v1.P[end] = pressure(v1.A[end], v1.A0[end], v1.beta[end], v1.Pext)
#   v2.P[1] = pressure(v2.A[1], v2.A0[1], v2.beta[1], v2.Pext)
#   v3.P[1] = pressure(v3.A[1], v3.A0[1], v3.beta[1], v3.Pext)
#
#   v1.c[end] = waveSpeed(v1.A[end], v1.gamma[end])
#   v2.c[1] = waveSpeed(v2.A[1], v2.gamma[1])
#   v3.c[1] = waveSpeed(v3.A[1], v3.gamma[1])
# end
#
# function calculateWstarBif(U :: Array, k :: Array)
#
#   W1 = U[1] + 4*k[1]*U[4]
#   W2 = U[2] - 4*k[2]*U[5]
#   W3 = U[3] - 4*k[3]*U[6]
#
#   return [W1, W2, W3]
# end
#
# function calculateFofUBif(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel,
#                            U :: Array,   k :: Array,   W :: Array)
#
#   f1 = U[1] + 4*k[1]*U[4] - W[1]
#
#   f2 = U[2] - 4*k[2]*U[5] - W[2]
#
#   f3 = U[3] - 4*k[3]*U[6] - W[3]
#
#   f4 = U[1]*(U[4]*U[4]*U[4]*U[4]) - U[2]*(U[5]*U[5]*U[5]*U[5]) - U[3]*(U[6]*U[6]*U[6]*U[6])
#
#   f5 = v1.beta[end]*(U[4]*U[4]/sqrt(v1.A0[end])  -1 ) -
#      ( v2.beta[1]*(U[5]*U[5]/sqrt(v2.A0[1]) - 1) )
#
#   f6 = v1.beta[end]*(U[4]*U[4]/sqrt(v1.A0[end])  -1 ) -
#      ( v3.beta[1]*(U[6]*U[6]/sqrt(v3.A0[1]) - 1) )
#
#   # f5 = v1.beta[end]*(U[4]*U[4]/sqrt(v1.A0[end])  -1 ) -
#   #    ( v2.beta[1]*(U[5]*U[5]/sqrt(v2.A0[1]) - 1) ) + 1060*0.5*(U[1]*U[1]-U[2]*U[2])
#   #
#   # f6 = v1.beta[end]*(U[4]*U[4]/sqrt(v1.A0[end])  -1 ) -
#   #    ( v3.beta[1]*(U[6]*U[6]/sqrt(v3.A0[1]) - 1) ) + 1060*0.5*(U[1]*U[1]-U[3]*U[3])
#
#   return [f1, f2, f3, f4, f5, f6]
# end
#
# function calculateJacobianBif(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel,
#                                U :: Array,   k :: Array)
#
#   J = eye(6)
#
#   J[1,4] =  4*k[1]
#   J[2,5] = -4*k[2]
#   J[3,6] = -4*k[3]
#
#   J[4,1] =  (U[4]*U[4]*U[4]*U[4])
#   J[4,2] = -(U[5]*U[5]*U[5]*U[5])
#   J[4,3] = -(U[6]*U[6]*U[6]*U[6])
#   J[4,4] =  4*U[1]*(U[4]*U[4]*U[4])
#   J[4,5] = -4*U[2]*(U[5]*U[5]*U[5])
#   J[4,6] = -4*U[3]*(U[6]*U[6]*U[6])
#
#   J[5,1] =  0.
#   # J[5,1] = 1060*U[1]
#
#   J[5,2] =  0.
#   # J[5,2] = -1060*U[2]
#
#   J[5,4] =  2*v1.beta[end]*U[4]/sqrt(v1.A0[end])
#   J[5,5] = -2*v2.beta[1]*U[5]/sqrt(v2.A0[1])
#
#   J[6,1] =  0.
#   # J[6,1] = 1060*U[1]
#
#   J[6,3] =  0.
#   # J[6,3] = -1060U[3]
#
#
#   J[6,4] =  2*v1.beta[end]*U[4]/sqrt(v1.A0[end])
#   J[6,6] = -2*v3.beta[1]*U[6]/sqrt(v3.A0[1])
#
#   return J
# end
function solveBifurcation(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel)

  U = [v1.u[end],
       v2.u[1],
       v3.u[1],
       sqrt(sqrt(v1.A[end])),
       sqrt(sqrt(v2.A[1])),
       sqrt(sqrt(v3.A[1]))]

  k1 = sqrt(1.5*v1.gamma[end])
  k2 = sqrt(1.5*v2.gamma[1])
  k3 = sqrt(1.5*v3.gamma[1])
  k = [k1, k2, k3]

  J = calculateJacobianBifurcation(v1, v2, v3, U, k)
  U = newtonRaphson(J, v1, v2, v3, U, k, calculateWstarBifurcation, calculateFBifurcation)

  updateBifurcation(U, v1, v2, v3)
end

function updateBifurcation(U :: Array, v1 :: Vessel, v2 :: Vessel, v3 :: Vessel)
    v1.u[end] = U[1]
    v2.u[1] = U[2]
    v3.u[1] = U[3]

    v1.A[end] = U[4]*U[4]*U[4]*U[4]
    v2.A[1] = U[5]*U[5]*U[5]*U[5]
    v3.A[1] = U[6]*U[6]*U[6]*U[6]

    v1.Q[end] = v1.u[end]*v1.A[end]
    v2.Q[1] = v2.u[1]*v2.A[1]
    v3.Q[1] = v3.u[1]*v3.A[1]

    v1.P[end] = pressure(v1.A[end], v1.A0[end], v1.beta[end], v1.Pext)
    v2.P[1] = pressure(v2.A[1], v2.A0[1], v2.beta[1], v2.Pext)
    v3.P[1] = pressure(v3.A[1], v3.A0[1], v3.beta[1], v3.Pext)

    v1.c[end] = waveSpeed(v1.A[end], v1.gamma[end])
    v2.c[1] = waveSpeed(v2.A[1], v2.gamma[1])
    v3.c[1] = waveSpeed(v3.A[1], v3.gamma[1])
end

function calculateWstarBifurcation(U :: Array, k :: Array)

  W1 = U[1] + 4*k[1]*U[4]
  W2 = U[2] - 4*k[2]*U[5]
  W3 = U[3] - 4*k[3]*U[6]

  return [W1, W2, W3]
end


function calculateJacobianBifurcation(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel,
                                      U :: Array, k :: Array)
  J = eye(6)

  J[1,4] =  4.0*k[1]
  J[2,5] = -4.0*k[2]
  J[3,6] = -4.0*k[3]

  J[4,1] =  U[4]*U[4]*U[4]*U[4]
  J[4,2] = -U[5]*U[5]*U[5]*U[5]
  J[4,3] = -U[6]*U[6]*U[6]*U[6]
  J[4,4] =  4.0*U[1]*(U[4]*U[4]*U[4])
  J[4,5] = -4.0*U[2]*(U[5]*U[5]*U[5])
  J[4,6] = -4.0*U[3]*(U[6]*U[6]*U[6])

  J[5,1] =  0.0
  # J[5,1] = 1060*U[1]

  J[5,2] =  0.0
  # J[5,2] = -1060*U[2]

  J[5,4] =  2.0*v1.beta[end]*U[4]*v1.s_inv_A0[end]
  J[5,5] = -2.0*v2.beta[1]*U[5]*v2.s_inv_A0[1]

  J[6,1] =  0.0
  # J[6,1] = 1060*U[1]

  J[6,3] =  0.0
  # J[6,3] = -1060U[3]

  J[6,4] =  2.0*v1.beta[end]*U[4]*v1.s_inv_A0[end]
  J[6,6] = -2.0*v3.beta[1]*U[6]*v3.s_inv_A0[1]

  return J
end

function calculateFBifurcation(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel,
                               U :: Array, k :: Array, W :: Array)

  f1 = U[1] + 4*k[1]*U[4] - W[1]

  f2 = U[2] - 4*k[2]*U[5] - W[2]

  f3 = U[3] - 4*k[3]*U[6] - W[3]

  f4 = U[1]*(U[4]*U[4]*U[4]*U[4]) - U[2]*(U[5]*U[5]*U[5]*U[5]) - U[3]*(U[6]*U[6]*U[6]*U[6])

  f5 = v1.beta[end]*(U[4]*U[4]*v1.s_inv_A0[end] - 1) -
      (v2.beta[1]*(U[5]*U[5]*v2.s_inv_A0[1] - 1))

  f6 = v1.beta[end]*(U[4]*U[4]*v1.s_inv_A0[end] - 1) -
      (v3.beta[1]*(U[6]*U[6]*v3.s_inv_A0[1] - 1))

  # f5 = v1.beta[end]*(U[4]*U[4]/sqrt(v1.A0[end])  -1 ) -
  #    ( v2.beta[1]*(U[5]*U[5]/sqrt(v2.A0[1]) - 1) ) + 1060*0.5*(U[1]*U[1]-U[2]*U[2])
  #
  # f6 = v1.beta[end]*(U[4]*U[4]/sqrt(v1.A0[end])  -1 ) -
  #    ( v3.beta[1]*(U[6]*U[6]/sqrt(v3.A0[1]) - 1) ) + 1060*0.5*(U[1]*U[1]-U[3]*U[3])

  return [f1, f2, f3, f4, f5, f6]
end
