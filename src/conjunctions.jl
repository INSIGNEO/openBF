function join_vessels!(b::Blood, v1::Vessel, v2::Vessel)

    U = SA[v1.u[end], v2.u[1], sqrt(sqrt(v1.A[end])), sqrt(sqrt(v2.A[1]))]

    k = SA[sqrt(1.5 * v1.gamma[end]), sqrt(1.5 * v2.gamma[1])]
    J = jacobian_conj(b, v1, v2, U, k)

    W = w_star_conj(U, k)
    F = f_conj(b, v1, v2, U, k, W)

    nr_toll_U = 1.e-5
    nr_toll_F = 1.e-5

    while true
        dU = J \ (-F)
        U += dU

        if any(isnan(dot(F, F)))
            println(F)
            @printf "error at conjunction with vessels %s and %s \n" v1.label v2.label
            break
        end

        u_ok = sum(abs.(dU) .<= nr_toll_U)
        f_ok = sum(abs.(F) .<= nr_toll_F)

        if u_ok == 4 || f_ok == 4
            break
        else
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
end


function w_star_conj(U::SVector, k::SVector)
    W1 = U[1] + 4 * k[1] * U[3]
    W2 = U[2] - 4 * k[2] * U[4]
    @SArray [W1, W2]
end


function f_conj(b::Blood, v1::Vessel, v2::Vessel, U::SVector, k::SVector, W::SVector)

    f1 = U[1] + 4 * k[1] * U[3] - W[1]

    f2 = U[2] - 4 * k[2] * U[4] - W[2]

    f3 = U[1] * (U[3] * U[3] * U[3] * U[3]) - U[2] * (U[4] * U[4] * U[4] * U[4])

    f4 =
        0.5 * b.rho * U[1] * U[1] + v1.beta[end] * (U[3] * U[3] / sqrt(v1.A0[end]) - 1) -
        (0.5 * b.rho * U[2] * U[2] + v2.beta[1] * (U[4] * U[4] / sqrt(v2.A0[1]) - 1))

    @SArray [f1, f2, f3, f4]
end


function jacobian_conj(b::Blood, v1::Vessel, v2::Vessel, U::SVector, k::SVector)

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

    SMatrix{4,4}(J)
end
