function join_vessels!(b::Blood, v1::Vessel, v2::Vessel)

    U = [v1.u[end], v2.u[1], sqrt(sqrt(v1.A[end])), sqrt(sqrt(v2.A[1]))]
    k = sqrt.(1.5 .* (v1.gamma[end], v2.gamma[1]))
    W = w_star_conj(U, k)
    J = jacobian_conj(b, v1, v2, U, k)
    F = f_conj(b, v1, v2, U, k, W)

    nr_toll_U = 1.e-5
    nr_toll_F = 1.e-5

    while true
        dU = J \ (-F)
        U_new = U + dU

        if any(isnan(dot(F, F)))
            println(F)
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
            W = w_star_conj(U, k)
            F = f_conj(b, v1, v2, U, k, W)
        end
    end

    v1.u[end] = U[1]
    v2.u[1] = U[2]

    v1.A[end] = U[3] * U[3] * U[3] * U[3]
    v1.Q[end] = v1.u[end] * v1.A[end]

    v2.A[1] = U[4]^4
    v2.Q[1] = v2.u[1] * v2.A[1]

    v1.P[end] = pressure(v1.A[end], v1.A0[end], v1.beta[end], v1.Pext)
    v2.P[1] = pressure(v2.A[1], v2.A0[1], v2.beta[1], v2.Pext)

    v1.c[end] = wave_speed(v1.A[end], v1.gamma[end])
    v2.c[1] = wave_speed(v2.A[1], v2.gamma[1])
end


function w_star_conj(U::Array, k::Tuple)

    W1 = U[1] + 4 * k[1] * U[3]
    W2 = U[2] - 4 * k[2] * U[4]

    [W1, W2]
end


function f_conj(b::Blood, v1::Vessel, v2::Vessel, U::Array, k::Tuple, W::Array)

    f1 = U[1] + 4 * k[1] * U[3] - W[1]

    f2 = U[2] - 4 * k[2] * U[4] - W[2]

    f3 = U[1] * (U[3] * U[3] * U[3] * U[3]) - U[2] * (U[4] * U[4] * U[4] * U[4])

    f4 =
        0.5 * b.rho * U[1] * U[1] + v1.beta[end] * (U[3] * U[3] / sqrt(v1.A0[end]) - 1) -
        (0.5 * b.rho * U[2] * U[2] + v2.beta[1] * (U[4] * U[4] / sqrt(v2.A0[1]) - 1))

    [f1, f2, f3, f4]
end


function jacobian_conj(b::Blood, v1::Vessel, v2::Vessel, U::Array, k::Tuple)

    J = zeros(4, 4) + I(4)

    J[1, 3] = 4 * k[1]
    J[2, 4] = -4 * k[2]

    J[3, 1] = U[3] * U[3] * U[3] * U[3]
    J[3, 2] = -U[4] * U[4] * U[4] * U[4]
    J[3, 3] = 4 * U[1] * (U[3] * U[3] * U[3]^3)
    J[3, 4] = -4 * U[2] * (U[4] * U[4] * U[4]^3)

    J[4, 3] = 2 * v1.beta[end] * U[3] / sqrt(v1.A0[end])
    J[4, 4] = -2 * v2.beta[1] * U[4] / sqrt(v2.A0[1])

    J[4, 1] = b.rho * U[1]
    J[4, 2] = -b.rho * U[2]

    J
end
