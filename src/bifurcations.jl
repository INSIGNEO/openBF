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


getUbif(v1::Vessel, v2::Vessel, v3::Vessel) = SVector{6,Float64}(v1.u[end], v2.u[1], v3.u[1],
        sqrt(sqrt(v1.A[end])), sqrt(sqrt(v2.A[1])), sqrt(sqrt(v3.A[1])))

function getJbif(v1::Vessel, v2::Vessel, v3::Vessel, U, k)
    J::Array{Float64, 2} = zeros(Float64, 6,6)
    
    J[1, 1] = 1.0
    J[2, 2] = 1.0
    J[3, 3] = 1.0
    J[4, 4] = 1.0

    J[1, 4] = 4k[1]
    J[2, 5] = -4k[2]
    J[3, 6] = -4k[3]

    J[4, 1] = (U[4] * U[4] * U[4] * U[4])
    J[4, 2] = -(U[5] * U[5] * U[5] * U[5])
    J[4, 3] = -(U[6] * U[6] * U[6] * U[6])
    J[4, 4] = 4U[1] * (U[4] * U[4] * U[4])
    J[4, 5] = -4U[2] * (U[5] * U[5] * U[5])
    J[4, 6] = -4U[3] * (U[6] * U[6] * U[6])

    J[5, 4] = 2v1.beta[end] * U[4] / sqrt(v1.A0[end])
    J[5, 5] = -2v2.beta[1] * U[5] / sqrt(v2.A0[1])

    J[6, 4] = 2v1.beta[end] * U[4] / sqrt(v1.A0[end])
    J[6, 6] = -2v3.beta[1] * U[6] / sqrt(v3.A0[1])

    SMatrix{6, 6, Float64, 36}(J)
end

function getF(v1::Vessel, v2::Vessel, v3::Vessel, U, k, W)
    SVector{6,Float64}(U[1] + 4k[1] * U[4] - W[1],
        U[2] - 4k[2] * U[5] - W[2],
        U[3] - 4k[3] * U[6] - W[3],
        U[1] * (U[4] * U[4] * U[4] * U[4]) - U[2] * (U[5] * U[5] * U[5] * U[5]) -
        U[3] * (U[6] * U[6] * U[6] * U[6]),
        v1.beta[end] * (U[4] * U[4] / sqrt(v1.A0[end]) - 1) -
        (v2.beta[1] * (U[5] * U[5] / sqrt(v2.A0[1]) - 1)),
        v1.beta[end] * (U[4] * U[4] / sqrt(v1.A0[end]) - 1) -
        (v3.beta[1] * (U[6] * U[6] / sqrt(v3.A0[1]) - 1)))
end

function NRbif(U, W, J, F, k, v1::Vessel, v2::Vessel, v3::Vessel)
    while norm(F)>1e-5
        U += J \ (-F)
        F = getF(v1, v2, v3, U, k, W)
        J = getJbif(v1, v2, v3, U, k)
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

function join_vessels!(v1::Vessel, v2::Vessel, v3::Vessel)
    k = (sqrt(1.5*v1.gamma[end]), sqrt(1.5*v2.gamma[1]), sqrt(1.5*v3.gamma[1]))
    U = getUbif(v1, v2, v3)
    W = (U[1] + 4k[1] * U[4], U[2] - 4k[2] * U[5], U[3] - 4k[3] * U[6])
    F = getF(v1, v2, v3, U, k, W)
    J = getJbif(v1, v2, v3, U, k)

    # solve
    U = NRbif(U, W, J, F, k, v1, v2, v3)

    updateBif!(v1, v2, v3, U)
end
