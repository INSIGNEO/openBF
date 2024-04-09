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

function solveAnastomosis(v1::Vessel, v2::Vessel, v3::Vessel)

    #Unknowns vector
    U = @SArray [
        v1.u[end],
        v2.u[end],
        v3.u[1],
        sqrt(sqrt(v1.A[end])),
        sqrt(sqrt(v2.A[end])),
        sqrt(sqrt(v3.A[1])),
    ]

    #Parameters vector
    k = @SArray [sqrt(1.5*v1.gamma[end]), sqrt(1.5*v2.gamma[end]), sqrt(1.5*v3.gamma[1])]
    W = calculateWstarAn(U, k)
    J = calculateJacobianAn(v1, v2, v3, U, k)
    F = calculateFofUAn(v1, v2, v3, U, k, W)

    #Newton-Raphson
    nr_toll_U = 1.e-5
    nr_toll_F = 1.e-5

    while true
        dU = J \ (-F)
        U_new = U + dU

        u_ok = 0
        f_ok = 0
        for i = 1:length(dU)
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
            W = calculateWstarBif(U, k)
            F = calculateFofUBif(v1, v2, v3, U, k, W)
        end
    end

    #Update vessel quantities
    v1.u[end] = U[1]
    v2.u[end] = U[2]
    v3.u[1] = U[3]

    v1.A[end] = U[4] * U[4] * U[4] * U[4]
    v2.A[end] = U[5] * U[5] * U[5] * U[5]
    v3.A[1] = U[6] * U[6] * U[6] * U[6]

    v1.Q[end] = v1.u[end] * v1.A[end]
    v2.Q[end] = v2.u[end] * v2.A[end]
    v3.Q[1] = v3.u[1] * v3.A[1]
end

function calculateWstarAn(U::SArray, k::SArray)

    W1 = U[1] + 4 * k[1] * U[4]
    W2 = U[2] + 4 * k[2] * U[5]
    W3 = U[3] - 4 * k[3] * U[6]

    return @SArray [W1, W2, W3]
end

function calculateFofUAn(v1::Vessel, v2::Vessel, v3::Vessel, U::SArray, k::SArray, W::SArray)

    f1 = U[1] + 4 * k[1] * U[4] - W[1]

    f2 = U[2] + 4 * k[2] * U[5] - W[2]

    f3 = U[3] - 4 * k[3] * U[6] - W[3]

    f4 =
        U[1] * (U[4] * U[4] * U[4] * U[4]) + U[2] * (U[5] * U[5] * U[5] * U[5]) -
        U[3] * (U[6] * U[6] * U[6] * U[6])

    f5 =
        v1.beta[end] * (U[4]^2 / sqrt(v1.A0[end]) - 1) -
        (v3.beta[1] * ((U[6])^2 / sqrt(v3.A0[1]) - 1))

    f6 =
        v2.beta[end] * (U[5]^2 / sqrt(v2.A0[end]) - 1) -
        (v3.beta[1] * ((U[6]^2) / sqrt(v3.A0[1]) - 1))

    return @SArray [f1, f2, f3, f4, f5, f6]
end

function calculateJacobianAn(v1::Vessel, v2::Vessel, v3::Vessel, U::SArray, k::SArray)
    J = @MArray zeros(6, 6) 
    J .+= I(6)

    J[1, 4] = 4 * k[1]
    J[2, 5] = 4 * k[2]
    J[3, 6] = -4 * k[3]

    J[4, 1] = U[4] * U[4] * U[4] * U[4]
    J[4, 2] = U[5] * U[5] * U[5] * U[5]
    J[4, 3] = -U[6] * U[6] * U[6] * U[6]
    J[4, 4] = 4 * U[1] * (U[4] * U[4] * U[4])
    J[4, 5] = 4 * U[2] * (U[5] * U[5] * U[5])
    J[4, 6] = -4 * U[3] * (U[6] * U[6] * U[6])

    J[5, 4] = 2 * v1.beta[end] * U[4] / sqrt(v1.A0[end])
    J[5, 5] = 0.0
    J[5, 6] = -2 * v3.beta[1] * U[6] / sqrt(v3.A0[1])

    J[6, 5] = 2 * v2.beta[end] * U[5] / sqrt(v2.A0[end])
    J[6, 6] = -2 * v3.beta[1] * U[6] / sqrt(v3.A0[1])

    J
end
