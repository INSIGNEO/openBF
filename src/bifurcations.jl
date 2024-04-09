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


function w!(W, U, k)
    W[1] = U[1] + 4 * k[1] * U[4]
    W[2] = U[2] - 4 * k[2] * U[5]
    W[3] = U[3] - 4 * k[3] * U[6]
end

function jac!(J, v1::Vessel, v2::Vessel, v3::Vessel, U, k)
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

    J
end

function f!(f, v1::Vessel, v2::Vessel, v3::Vessel, U, k, W)
    f[1] = U[1] + 4 * k[1] * U[4] - W[1]
    f[2] = U[2] - 4 * k[2] * U[5] - W[2]
    f[3] = U[3] - 4 * k[3] * U[6] - W[3]
    f[4] =
        U[1] * (U[4] * U[4] * U[4] * U[4]) - U[2] * (U[5] * U[5] * U[5] * U[5]) -
        U[3] * (U[6] * U[6] * U[6] * U[6])
    f[5] =
        v1.beta[end] * (U[4] * U[4] / sqrt(v1.A0[end]) - 1) -
        (v2.beta[1] * (U[5] * U[5] / sqrt(v2.A0[1]) - 1))
    f[6] =
        v1.beta[end] * (U[4] * U[4] / sqrt(v1.A0[end]) - 1) -
        (v3.beta[1] * (U[6] * U[6] / sqrt(v3.A0[1]) - 1))
end

function u!(U, v1::Vessel, v2::Vessel, v3::Vessel)
    U[1] = v1.u[end]
    U[2] = v2.u[1]
    U[3] = v3.u[1]
    U[4] = sqrt(sqrt(v1.A[end]))
    U[5] = sqrt(sqrt(v2.A[1]))
    U[6] = sqrt(sqrt(v3.A[1]))
end

function join_vessels!(n::Network, v1::Vessel, v2::Vessel, v3::Vessel)
    u!(n.bifU, v1, v2, v3)

    k = (sqrt(1.5*v1.gamma[end]), sqrt(1.5*v2.gamma[1]), sqrt(1.5*v3.gamma[1]))
    
    w!(n.bifW, n.bifU, k)
    jac!(n.bifJ, v1, v2, v3, n.bifU, k)
    
    f!(n.bifF, v1, v2, v3, n.bifU, k, n.bifW)

    # solve
    while norm(n.bifF) > 1e-5
        n.bifU .+= n.bifJ \ (-n.bifF)
        w!(n.bifW, n.bifU, k)
        f!(n.bifF, v1, v2, v3, n.bifU, k, n.bifW)
    end

    # update
    v1.u[end] = n.bifU[1]
    v2.u[1] = n.bifU[2]
    v3.u[1] = n.bifU[3]

    v1.A[end] = n.bifU[4] * n.bifU[4] * n.bifU[4] * n.bifU[4]
    v2.A[1] = n.bifU[5] * n.bifU[5] * n.bifU[5] * n.bifU[5]
    v3.A[1] = n.bifU[6] * n.bifU[6] * n.bifU[6] * n.bifU[6]

    v1.Q[end] = v1.u[end] * v1.A[end]
    v2.Q[1] = v2.u[1] * v2.A[1]
    v3.Q[1] = v3.u[1] * v3.A[1]
end
