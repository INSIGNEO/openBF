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

getUconj(v1::Vessel, v2::Vessel) = SVector{4,Float64}(v1.u[end], v2.u[1],
    sqrt(sqrt(v1.A[end])), sqrt(sqrt(v2.A[1])))

getWconj(U, k) = SVector{2,Float64}(U[1] + 4k[1] * U[3], U[2] - 4k[2] * U[4])

getFconj(v1::Vessel, v2::Vessel, U, k, W, ρ) = SVector{4,Float64}(U[1] + 4k[1] * U[3] - W[1],
        U[2] - 4k[2] * U[4] - W[2],
        U[1] * (U[3] * U[3] * U[3] * U[3]) - U[2] * (U[4] * U[4] * U[4] * U[4]),
        0.5 * ρ * U[1] * U[1] + v1.beta[end] * (U[3] * U[3] / sqrt(v1.A0[end]) - 1) -
        (0.5 * ρ * U[2] * U[2] + v2.beta[1] * (U[4] * U[4] / sqrt(v2.A0[1]) - 1)))

function getJconj(v1::Vessel, v2::Vessel, U, k, ρ)
    J::Array{Float64, 2} = zeros(Float64, 4, 4)
    
    J[1, 1] = 1.0
    J[2, 2] = 1.0

    J[1, 3] = 4k[1]
    J[2, 4] = -4k[2]

    J[3, 1] = U[3] * U[3] * U[3] * U[3]
    J[3, 2] = -U[4] * U[4] * U[4] * U[4]
    J[3, 3] = 4U[1] * U[3] * U[3] * U[3]
    J[3, 4] = -4U[2] * U[4] * U[4] * U[4]

    J[4, 3] = 2v1.beta[end] * U[3] / sqrt(v1.A0[end])
    J[4, 4] = -2v2.beta[1] * U[4] / sqrt(v2.A0[1])

    J[4, 1] = ρ * U[1]
    J[4, 2] = -ρ * U[2]

    SMatrix{4, 4, Float64, 16}(J)
end

function NRconj(U, W, J, F, k, v1::Vessel, v2::Vessel, ρ)
    while norm(F)>1e-5
        U += J \ (-F)
        W = getWconj(U, k)
        F = getFconj(v1, v2, U, k, W, ρ)
        J = getJconj(v1, v2, U, k, ρ)
    end
    U
end

function updateConj!(v1::Vessel, v2::Vessel, U)
    v1.u[end] = U[1]
    v2.u[1] = U[2]

    v1.A[end] = U[3] * U[3] * U[3] * U[3]
    v1.Q[end] = v1.u[end] * v1.A[end]

    v2.A[1] = U[4] * U[4] * U[4] * U[4]
    v2.Q[1] = v2.u[1] * v2.A[1]
end

function join_vessels!(v1::Vessel, v2::Vessel, ρ::Float64)
    k = (sqrt(1.5*v1.gamma[end]), sqrt(1.5*v2.gamma[1]))
    U = getUconj(v1, v2)
    W = getWconj(U, k)
    F = getFconj(v1, v2, U, k, W, ρ)
    J = getJconj(v1, v2, U, k, ρ)

    # solve
    U = NRconj(U, W, J, F, k, v1, v2, ρ)

    updateConj!(v1, v2, U)
end
