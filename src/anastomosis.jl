
function join_anastomosis!(v1::Vessel, v2::Vessel, v3::Vessel)

    #Unknowns vector
    U = SA[
        v1.u[end],
        v2.u[end],
        v3.u[1],
        sqrt(sqrt(v1.A[end])),
        sqrt(sqrt(v2.A[end])),
        sqrt(sqrt(v3.A[1])),
    ]

    #Parameters vector
    k = SA[sqrt(1.5 * v1.gamma), sqrt(1.5 * v2.gamma), sqrt(1.5 * v3.gamma)]
    J = jacobian_anastomosis(v1, v2, v3, U, k)

    W = w_star_anastomosis(U, k)
    F = f_anastomosis(v1, v2, v3, U, k, W)

    #Newton-Raphson
    nr_toll_U = 1.e-5
    nr_toll_F = 1.e-5

    while true
        dU = J \ (-F)
        U += dU

        if any(isnan(dot(F, F)))
            println(F)
            @printf "error at anastomosis with vessels %s, %s, and %s \n" v1.label v2.label v3.label
            break
        end

        u_ok = sum(abs.(dU) .<= nr_toll_U)
        f_ok = sum(abs.(F) .<= nr_toll_F)

        if u_ok == 6 || f_ok == 6
            break
        else
            W = w_star_anastomosis(U, k)
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

    v1.P[end] = pressure(v1.A[end], v1.A0[end], v1.beta[end], v1.Pext)
    v2.P[end] = pressure(v2.A[end], v2.A0[end], v2.beta[end], v2.Pext)
    v3.P[1] = pressure(v3.A[1], v3.A0[1], v3.beta[1], v3.Pext)

    v1.c[end] = wave_speed(v1.A[end], v1.gamma[end])
    v2.c[end] = wave_speed(v2.A[end], v2.gamma[end])
    v3.c[1] = wave_speed(v3.A[1], v3.gamma[1])

end

function w_star_anastomosis(U::SVector, k::SVector)
    W1 = U[1] + 4k[1] * U[4]
    W2 = U[2] + 4k[2] * U[5]
    W3 = U[3] - 4k[3] * U[6]

    @SArray [W1, W2, W3]
end

function f_anastomosis(
    v1::Vessel,
    v2::Vessel,
    v3::Vessel,
    U::SVector,
    k::SVector,
    W::SVector,
)

    f1 = U[1] + 4k[1] * U[4] - W[1]
    f2 = U[2] + 4k[2] * U[5] - W[2]
    f3 = U[3] - 4k[3] * U[6] - W[3]
    f4 =
        U[1] * (U[4] * U[4] * U[4] * U[4]) + U[2] * (U[5] * U[5] * U[5] * U[5]) -
        U[3] * (U[6] * U[6] * U[6] * U[6])
    f5 =
        v1.beta[end] * (U[4]^2 / sqrt(v1.A0[end]) - 1) -
        (v3.beta[1] * ((U[6])^2 / sqrt(v3.A0[1]) - 1))
    f6 =
        v2.beta[end] * (U[5]^2 / sqrt(v2.A0[end]) - 1) -
        (v3.beta[1] * ((U[6]^2) / sqrt(v3.A0[1]) - 1))

    @SArray [f1, f2, f3, f4, f5, f6]
end

function jacobian_anastomosis(v1::Vessel, v2::Vessel, v3::Vessel, U::Array, k::Array)
    J = zeros(6, 6) + I(6)

    J[1, 4] = 4k[1]
    J[2, 5] = 4k[2]
    J[3, 6] = -4k[3]

    J[4, 1] = U[4] * U[4] * U[4] * U[4]
    J[4, 2] = U[5] * U[5] * U[5] * U[5]
    J[4, 3] = -U[6] * U[6] * U[6] * U[6]
    J[4, 4] = 4U[1] * (U[4] * U[4] * U[4])
    J[4, 5] = 4U[2] * (U[5] * U[5] * U[5])
    J[4, 6] = -4U[3] * (U[6] * U[6] * U[6])

    J[5, 4] = 2v1.beta[end] * U[4] / sqrt(v1.A0[end])
    J[5, 5] = 0.0
    J[5, 6] = -2v3.beta[1] * U[6] / sqrt(v3.A0[1])

    J[6, 5] = 2v2.beta[end] * U[5] / sqrt(v2.A0[end])
    J[6, 6] = -2v3.beta[1] * U[6] / sqrt(v3.A0[1])

    SMatrix{6,6}(J)
end
