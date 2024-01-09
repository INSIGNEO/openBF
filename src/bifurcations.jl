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

function join_vessels!(v1::Vessel, v2::Vessel, v3::Vessel)
    #Unknowns vector
    U = [
        v1.u[end],
        v2.u[1],
        v3.u[1],
        sqrt(sqrt(v1.A[end])),
        sqrt(sqrt(v2.A[1])),
        sqrt(sqrt(v3.A[1])),
    ]

    #Parameters vector
    k = sqrt.(1.5 .* (v1.gamma[end], v2.gamma[1], v3.gamma[1]))
    W = w_star_bifurcation(U, k)
    J = jacobian_bifurcation(v1, v2, v3, U, k)
    F = f_bifurcation(v1, v2, v3, U, k, W)

    #Newton-Raphson
    nr_toll_U = 1e-5
    nr_toll_F = 1e-5

    while true
        dU = J \ (-F)
        # U_new = U + 0.01*dU
        U_new = U + dU

        if any(isnan(dot(F, F)))
            println(F)
            @printf "error at bifurcation with vessels %s, %s, and %s \n" v1.label v2.label v3.label
            break
        end

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
            W = w_star_bifurcation(U, k)
            F = f_bifurcation(v1, v2, v3, U, k, W)
        end
    end

    #Update vessel quantities
    v1.u[end] = U[1]
    v2.u[1] = U[2]
    v3.u[1] = U[3]

    v1.A[end] = U[4] * U[4] * U[4] * U[4]
    v2.A[1] = U[5] * U[5] * U[5] * U[5]
    v3.A[1] = U[6] * U[6] * U[6] * U[6]

    v1.Q[end] = v1.u[end] * v1.A[end]
    v2.Q[1] = v2.u[1] * v2.A[1]
    v3.Q[1] = v3.u[1] * v3.A[1]

    v1.P[end] = pressure(v1.A[end], v1.A0[end], v1.beta[end], v1.Pext)
    v2.P[1] = pressure(v2.A[1], v2.A0[1], v2.beta[1], v2.Pext)
    v3.P[1] = pressure(v3.A[1], v3.A0[1], v3.beta[1], v3.Pext)

    v1.c[end] = wave_speed(v1.A[end], v1.gamma[end])
    v2.c[1] = wave_speed(v2.A[1], v2.gamma[1])
    v3.c[1] = wave_speed(v3.A[1], v3.gamma[1])
end


function w_star_bifurcation(U::Array, k::Tuple)
    W1 = U[1] + 4 * k[1] * U[4]
    W2 = U[2] - 4 * k[2] * U[5]
    W3 = U[3] - 4 * k[3] * U[6]

    [W1, W2, W3]
end


function f_bifurcation(v1::Vessel, v2::Vessel, v3::Vessel, U::Array, k::Tuple, W::Array)
    f1 = U[1] + 4 * k[1] * U[4] - W[1]
    f2 = U[2] - 4 * k[2] * U[5] - W[2]
    f3 = U[3] - 4 * k[3] * U[6] - W[3]
    f4 =
        U[1] * (U[4] * U[4] * U[4] * U[4]) - U[2] * (U[5] * U[5] * U[5] * U[5]) -
        U[3] * (U[6] * U[6] * U[6] * U[6])
    f5 =
        v1.beta[end] * (U[4] * U[4] / sqrt(v1.A0[end]) - 1) -
        (v2.beta[1] * (U[5] * U[5] / sqrt(v2.A0[1]) - 1))
    f6 =
        v1.beta[end] * (U[4] * U[4] / sqrt(v1.A0[end]) - 1) -
        (v3.beta[1] * (U[6] * U[6] / sqrt(v3.A0[1]) - 1))
    [f1, f2, f3, f4, f5, f6]
end


function jacobian_bifurcation(v1::Vessel, v2::Vessel, v3::Vessel, U::Array, k::Tuple)
    J = zeros(6, 6) + I(6)

    J[1, 4] = 4 * k[1]
    J[2, 5] = -4 * k[2]
    J[3, 6] = -4 * k[3]

    J[4, 1] = (U[4] * U[4] * U[4] * U[4])
    J[4, 2] = -(U[5] * U[5] * U[5] * U[5])
    J[4, 3] = -(U[6] * U[6] * U[6] * U[6])
    J[4, 4] = 4 * U[1] * (U[4] * U[4] * U[4])
    J[4, 5] = -4 * U[2] * (U[5] * U[5] * U[5])
    J[4, 6] = -4 * U[3] * (U[6] * U[6] * U[6])

    J[5, 1] = 0.0
    J[5, 2] = 0.0
    J[5, 4] = 2 * v1.beta[end] * U[4] / sqrt(v1.A0[end])
    J[5, 5] = -2 * v2.beta[1] * U[5] / sqrt(v2.A0[1])

    J[6, 1] = 0.0
    J[6, 3] = 0.0
    J[6, 4] = 2 * v1.beta[end] * U[4] / sqrt(v1.A0[end])
    J[6, 6] = -2 * v3.beta[1] * U[6] / sqrt(v3.A0[1])

    J
end
