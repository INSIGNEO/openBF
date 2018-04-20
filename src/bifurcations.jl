# #=
# Copyright 2018 INSIGNEO Institute for in silico Medicine
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


"""
    solveBifurcation(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel)

Solve the non-linear system at the bifurcation node between three vessels.
"""
function solveBifurcation(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel)
    U0 = @SArray [v1.u[end],
       v2.u[1],
       v3.u[1],
       sqrt(sqrt(v1.A[end])),
       sqrt(sqrt(v2.A[1])),
       sqrt(sqrt(v3.A[1]))]

    k1 = v1.s_15_gamma[end]
    k2 = v2.s_15_gamma[1]
    k3 = v3.s_15_gamma[1]
    k = @SArray [k1, k2, k3]

<<<<<<< HEAD
    J = calculateJacobianBifurcation(v1, v2, v3, U0, k)
    U = newtonRaphson([v1, v2, v3], J, U0, k, calculateWstarBifurcation, calculateFBifurcation)

    updateBifurcation(U, v1, v2, v3)
end

=======
# *function* __`calculateWstarBif`__ $\rightarrow$
# `W::Array{Float, 1}`
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------- -------------------------------------------------------------
# `U`            `::Array` junction unknown vector.
#
# `k`            `::Array` junction parameters vector.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# Riemann invariants are computed starting from unknown and parameters
# vectors as
# $$
#   W_1 = u_1 + 4c_1, \quad W_2 = u_2 - 4c_2, \quad W_3 = u_3 - 4c_3.
# $$
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ----------- ----------------------------------------------------------------
# `W`         `::Array` containing outgoing characteristics from the three
#             vessels.
# ----------------------------------------------------------------------------
# <a name="calculateWstarBif"></a>
function calculateWstarBif(U :: Array{Float64, 1}, k :: Array{Float64, 1})

  W1 = U[1] + 4*k[1]*U[4]
  W2 = U[2] - 4*k[2]*U[5]
  W3 = U[3] - 4*k[3]*U[6]

  return [W1, W2, W3]
end

# *function* __`calculateFofUBif`__ $\rightarrow$ `F::Array`
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------- -------------------------------------------------------------
# `b`            `::Blood` data structure.
#
# `v1`           `::Vessel` parent vessel data structure.
#
# `v2`           `::Vessel` first daughter vessel data structure.
#
# `v3`           `::Vessel` second daughter vessel data structure.
#
# `U`            `::Array` junction unknown vector.
#
# `k`            `::Array` junction parameters vector.
#
# `W`            `::Array` outgoing characteristics array.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# $F$ array is computed by imposing the conservation of mass and static
# pressure at bifurcation node. $F$ reads
# $$
#   F = \left\{f_i \right\} =
#   \begin{cases}
#     U_1 + 4k_1U_4 - W_1^* = 0, \\
#     U_2 - 4k_2U_5 - W_2^* = 0, \\
#     U_3 - 4k_3U_6 - W_3^* = 0, \\
#     U_1U_4^4 - U_2U_5^4 - U_3U_6^4 = 0, \\
#     \beta_1 \left(\tfrac{U_4^2}{A_{01}^{1/2}} -1 \right) -
#     \beta_2 \left(\tfrac{U_5^2}{A_{02}^{1/2}} -1 \right) = 0, \\
#     \beta_1 \left(\tfrac{U_4^2}{A_{01}^{1/2}} -1 \right) -
#     \beta_3 \left(\tfrac{U_6^2}{A_{03}^{1/2}} -1 \right) = 0, \\
#   \end{cases}
# $$
# ----------------------------------------------------------------------------
# Total pressure conservation can be imposed by replacing `f5` and `f6` with
#
#     f5 = v1.beta*( U[4]^2 - sqrt(v1.A0) ) -
#        ( v2.beta*( U[5]^2 - sqrt(v2.A0) )
#
#     f6 = v1.beta*( U[4]^2 - sqrt(v1.A0) ) -
#        ( v3.beta*( U[6]^2 - sqrt(v3.A0) )
#
# respectively.
#
# ----------------------------------------------------------------------------
# Returns:
# ----------- ----------------------------------------------------------------
# `F`         `::Array` Newton's method relations.
# ----------------------------------------------------------------------------
# <a name="calculateFofUBif"></a>
function calculateFofUBif(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel,
                           U :: Array{Float64, 1},   k :: Array{Float64, 1},   W :: Array{Float64, 1})
>>>>>>> master

"""
    calculateJacobianBifurcation(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel, U, k)

Return the Jacobian for bifurcation equations.
"""
function calculateJacobianBifurcation(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel, U, k)
    U43 = U[4]*U[4]*U[4]
    U53 = U[5]*U[5]*U[5]
    U63 = U[6]*U[6]*U[6]

    J14 = 4.0*k[1]
    J25 = -4.0*k[2]
    J36 = -4.0*k[3]

    J41 = U[4]*U43
    J42 = -U[5]*U53
    J43 = -U[6]*U63
    J44 = 4.0*U[1]*U43
    J45 = -4.0*U[2]*U53
    J46 = -4.0*U[3]*U63

    J54 = 2.0*v1.beta[end]*U[4]*v1.s_inv_A0[end]
    J55 = -2.0*v2.beta[1]*U[5]*v2.s_inv_A0[1]

    J64 = 2.0*v1.beta[end]*U[4]*v1.s_inv_A0[end]
    J66 = -2.0*v3.beta[1]*U[6]*v3.s_inv_A0[1]

    return @SArray [1.0 0.0 0.0 J14 0.0 0.0;
                    0.0 1.0 0.0 0.0 J25 0.0;
                    0.0 0.0 1.0 0.0 0.0 J36;
                    J41 J42 J43 J44 J45 J46;
                    0.0 0.0 0.0 J54 J55 0.0;
                    0.0 0.0 0.0 J64 0.0 J66]
end


"""
    calculateWstarBifurcation(U, k)

Return the Riemann invariants at the bifurcation node.
"""
function calculateWstarBifurcation(U, k)
    W1 = U[1] + 4.0*k[1]*U[4]
    W2 = U[2] - 4.0*k[2]*U[5]
    W3 = U[3] - 4.0*k[3]*U[6]

    return @SArray [W1, W2, W3]
end

<<<<<<< HEAD
=======
# *function* __`calculateJacobianBif`__ $\rightarrow$
# `W::Array{Float, 1}`
#
# ----------------------------------------------------------------------------
# Parameters:
# -------------- -------------------------------------------------------------
# `b`            `::Blood` data structure.
#
# `v1`           `::Vessel` parent vessel data structure.
#
# `v2`           `::Vessel` daughter vessel data structure.
#
# `U`            `::Array` junction unknown vector.
#
# `k`            `::Array` junction parameters vector.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# The Jacobian is computed as
# $$
#   J = \left[ \begin{array}{cccccc}
#         1 & 0 & 0 & 4k_1 & 0 & 0 \\
#         0 & 1 & 0 & 0 & -4k_2 & 0 \\
#         0 & 0 & 1 & 0 & 0 & -4k_3 \\
#         U_4^4 & -U_5^4 & -U_6^4 & 4U_1 U_4^3 & -4 U_2 U_5^3 & -4U_3 U_6^3 \\
#         0 & 0 & 0 & 2\beta_1 U_4/A_{01}^{1/2} &
#                    -2\beta_2 U_5/A_{02}^{1/2} & 0 \\
#         0 & 0 & 0 & 2\beta_1 U_4/A_{01}^{1/2} &
#                    -2\beta_2 U_5/A_{02}^{1/2} & 0 \\
#       \end{array} \right].
# $$
# ----------------------------------------------------------------------------
#
# The conservation of total pressure is imposed by setting the following
# elements different than zero:
#
#     J[5,1] =  rho*U[1]
#     J[5,2] = -rho*U[2]
#     J[5,4] =  2*v1.beta*U[4]
#     J[5,5] = -2*v2.beta*U[5]
#
#     J[6,1] =  rho*U[1]
#     J[6,3] = -rho*U[3]
#     J[6,4] =  2*v1.beta*U[4]
#     J[6,6] = -2*v3.beta*U[6]
#
# ----------------------------------------------------------------------------
# Returns:
# ----------- ----------------------------------------------------------------
# `J`         `::Array{Float, 2}` Jacobian matrix.
# ----------------------------------------------------------------------------
# <a name="calculateJacobianBif"></a>
function calculateJacobianBif(v1 :: Vessel, v2 :: Vessel, v3 :: Vessel,
                               U :: Array{Float64, 1},   k :: Array{Float64, 1})
>>>>>>> master

"""
    calculateFBifurcation(vessels :: Array{Vessel,1}, U, k, W)
"""
function calculateFBifurcation(vessels :: Array{Vessel,1}, U, k, W)
    v1 = vessels[1]
    v2 = vessels[2]
    v3 = vessels[3]

    U42 = U[4]*U[4]
    U52 = U[5]*U[5]
    U62 = U[6]*U[6]

    f1 = U[1] + 4*k[1]*U[4] - W[1]

<<<<<<< HEAD
    f2 = U[2] - 4*k[2]*U[5] - W[2]
=======
  three_times_U4 = U[4]*U[4]*U[4]
  three_times_U5 = U[4]*U[4]*U[4]
  three_times_U6 = U[4]*U[4]*U[4]

  J[4,1] =   three_times_U4*U[4]
  J[4,2] = -(three_times_U5*U[5])
  J[4,3] = -(three_times_U6*U[6])
  J[4,4] =  4*U[1]*(three_times_U4)
  J[4,5] = -4*U[2]*(three_times_U5)
  J[4,6] = -4*U[3]*(three_times_U6)
>>>>>>> master

    f3 = U[3] - 4*k[3]*U[6] - W[3]

    f4 = U[1]*(U42*U42) - U[2]*(U52*U52) - U[3]*(U62*U62)

    f5 = v1.beta[end]*(U42*v1.s_inv_A0[end] - 1.0) -
        (v2.beta[1]*(U52*v2.s_inv_A0[1] - 1.0))

    f6 = v1.beta[end]*(U42*v1.s_inv_A0[end] - 1.0) -
        (v3.beta[1]*(U62*v3.s_inv_A0[1] - 1.0))

    return @SArray [f1, f2, f3, f4, f5, f6]
end


<<<<<<< HEAD
"""
    updateBifurcation(U, v1 :: Vessel, v2 :: Vessel, v3 :: Vessel)
=======
  J[5,4] =  2*v1.beta[end]*U[4]*v1.s_inv_A0[end]
  J[5,5] = -2*v2.beta[ 1 ]*U[5]*v2.s_inv_A0[1]
>>>>>>> master

Update the values at the bifurcation node for the three vessels.
"""
function updateBifurcation(U, v1 :: Vessel, v2 :: Vessel, v3 :: Vessel)
    v1.u[end] = U[1]
    v2.u[1] = U[2]
    v3.u[1] = U[3]

    v1.A[end] = U[4]*U[4]*U[4]*U[4]
    v2.A[1] = U[5]*U[5]*U[5]*U[5]
    v3.A[1] = U[6]*U[6]*U[6]*U[6]

<<<<<<< HEAD
    v1.Q[end] = v1.u[end]*v1.A[end]
    v2.Q[1] = v2.u[1]*v2.A[1]
    v3.Q[1] = v3.u[1]*v3.A[1]

    v1.P[end] = pressure(v1.A[end], v1.A0[end], v1.beta[end], v1.Pext)
    v2.P[1] = pressure(v2.A[1], v2.A0[1], v2.beta[1], v2.Pext)
    v3.P[1] = pressure(v3.A[1], v3.A0[1], v3.beta[1], v3.Pext)
=======
  J[6,4] =  J[5,4]
  J[6,6] = -2*v3.beta[ 1 ]*U[6]*v3.s_inv_A0[1]
>>>>>>> master

    v1.c[end] = waveSpeed(v1.A[end], v1.gamma[end])
    v2.c[1] = waveSpeed(v2.A[1], v2.gamma[1])
    v3.c[1] = waveSpeed(v3.A[1], v3.gamma[1])
end
