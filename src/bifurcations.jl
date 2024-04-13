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


function getW(U, k)
    SVector{3,Float64}(U[1] + 4 * k[1] * U[4], U[2] - 4 * k[2] * U[5], U[3] - 4 * k[3] * U[6])
end

function getJ(v1::Vessel, v2::Vessel, v3::Vessel, U, k)
    J = zeros(6,6)
    
    J[1, 1] = 1.0
    J[2, 2] = 1.0
    J[3, 3] = 1.0
    J[4, 4] = 1.0

    J[1, 4] = 4 * k[1]
    J[2, 5] = -4 * k[2]
    J[3, 6] = -4 * k[3]

    J[4, 1] = (U[4] * U[4] * U[4] * U[4])
    J[4, 2] = -(U[5] * U[5] * U[5] * U[5])
    J[4, 3] = -(U[6] * U[6] * U[6] * U[6])
    J[4, 4] = 4 * U[1] * (U[4] * U[4] * U[4])
    J[4, 5] = -4 * U[2] * (U[5] * U[5] * U[5])
    J[4, 6] = -4 * U[3] * (U[6] * U[6] * U[6])

    J[5, 4] = 2 * v1.beta[end] * U[4] / sqrt(v1.A0[end])
    J[5, 5] = -2 * v2.beta[1] * U[5] / sqrt(v2.A0[1])

    J[6, 4] = 2 * v1.beta[end] * U[4] / sqrt(v1.A0[end])
    J[6, 6] = -2 * v3.beta[1] * U[6] / sqrt(v3.A0[1])

    SMatrix{6, 6}(J)
end

function getF(v1::Vessel, v2::Vessel, v3::Vessel, U, k, W)
    SVector{6,Float64}(U[1] + 4 * k[1] * U[4] - W[1],
        U[2] - 4 * k[2] * U[5] - W[2],
        U[3] - 4 * k[3] * U[6] - W[3],
        U[1] * (U[4] * U[4] * U[4] * U[4]) - U[2] * (U[5] * U[5] * U[5] * U[5]) -
        U[3] * (U[6] * U[6] * U[6] * U[6]),
        v1.beta[end] * (U[4] * U[4] / sqrt(v1.A0[end]) - 1) -
        (v2.beta[1] * (U[5] * U[5] / sqrt(v2.A0[1]) - 1)),
        v1.beta[end] * (U[4] * U[4] / sqrt(v1.A0[end]) - 1) -
        (v3.beta[1] * (U[6] * U[6] / sqrt(v3.A0[1]) - 1)))
end

function getU(v1::Vessel, v2::Vessel, v3::Vessel)
    SVector{6,Float64}(v1.u[end], v2.u[1], v3.u[1],
        sqrt(sqrt(v1.A[end])), sqrt(sqrt(v2.A[1])), sqrt(sqrt(v3.A[1])))
end

function getU(U::SVector{6, Float64}, dU)
    SVector{6,Float64}(U[1]+dU[1], U[2]+dU[2], U[3]+dU[3], U[4]+dU[4], U[5]+dU[5], U[6]+dU[6])
end


function NR(U, W, J, F, k, v1::Vessel, v2::Vessel, v3::Vessel)
    while norm(F)>1e-5
        dU = J \ (-F)
        U = getU(U, dU)
        W = getW(U, k)
        F = getF(v1, v2, v3, U, k, W)
    end
    U
end

function updateBif!(v1::Vessel, v2::Vessel, v3::Vessel, U)
    v1.u[end] = U[1]
    v2.u[1] = U[2]
    v3.u[1] = U[3]

    v1.A[end] = U[4] * U[4] * U[4] * U[4]
    v2.A[1] = U[5] * U[5] * U[5] * U[5]
    v3.A[1] = U[6] * U[6] * U[6] * U[6]

    v1.Q[end] = v1.u[end] * v1.A[end]
    v2.Q[1] = v2.u[1] * v2.A[1]
    v3.Q[1] = v3.u[1] * v3.A[1]
end

function join_vessels!(n::Network, v1::Vessel, v2::Vessel, v3::Vessel)
    k = (sqrt(1.5*v1.gamma[end]), sqrt(1.5*v2.gamma[1]), sqrt(1.5*v3.gamma[1]))
    U = getU(v1, v2, v3)
    W = getW(U, k)
    J = getJ(v1, v2, v3, U, k)
    F = getF(v1, v2, v3, U, k, W)

    # solve
    U = NR(U, W, J, F, k, v1, v2, v3)

    updateBif!(v1, v2, v3, U)
end
