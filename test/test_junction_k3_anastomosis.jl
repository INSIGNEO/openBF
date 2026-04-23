# Validate solve_junction! at k=3 against the existing solveAnastomosis solver.
# Anastomosis uses static pressure (use_total_pressure=false).

@testset "Junction solve k=3 (anastomosis)" begin
    blood = openBF.Blood(Dict("rho" => 1060.0, "mu" => 0.004))
    cfg_p1 = Dict{Any,Any}(
        "label" => "parent1", "sn" => 1, "tn" => 3,
        "L" => 0.1, "E" => 400_000.0, "R0" => 4.0e-3, "h0" => 4e-4,
        "gamma_profile" => 9,
    )
    cfg_p2 = Dict{Any,Any}(
        "label" => "parent2", "sn" => 2, "tn" => 3,
        "L" => 0.1, "E" => 400_000.0, "R0" => 3.5e-3, "h0" => 3.5e-4,
        "gamma_profile" => 9,
    )
    cfg_d = Dict{Any,Any}(
        "label" => "daughter", "sn" => 3, "tn" => 4,
        "L" => 0.1, "E" => 500_000.0, "R0" => 5.5e-3, "h0" => 5.5e-4,
        "gamma_profile" => 9,
        "R1" => 5e7, "R2" => 2e9, "Cc" => 5e-10,
        "inlet_impedance_matching" => false,
    )

    v1 = openBF.Vessel(cfg_p1, blood, 100, ["P", "Q", "A", "u"])
    v2 = openBF.Vessel(cfg_p2, blood, 100, ["P", "Q", "A", "u"])
    v3 = openBF.Vessel(cfg_d,  blood, 100, ["P", "Q", "A", "u"])

    # Perturb boundary faces to a non-trivial state.
    v1.A[end] *= 1.02
    v1.Q[end]  = v1.A[end] * 0.20
    v1.u[end]  = v1.Q[end] / v1.A[end]

    v2.A[end] *= 1.01
    v2.Q[end]  = v2.A[end] * 0.18
    v2.u[end]  = v2.Q[end] / v2.A[end]

    v3.A[1]   *= 0.99
    v3.Q[1]    = v3.A[1] * 0.35
    v3.u[1]    = v3.Q[1] / v3.A[1]

    v1_leg = deepcopy(v1)
    v2_leg = deepcopy(v2)
    v3_leg = deepcopy(v3)

    # Legacy anastomosis solve.
    openBF.solveAnastomosis(v1_leg, v2_leg, v3_leg)

    # Generic junction solve (static pressure, matching anastomosis).
    jc = Junction(3, [1, 2, 3], [:outlet, :outlet, :inlet])
    solve_junction!(jc, [v1, v2, v3], blood)

    rtol = 1e-10

    @test isapprox(v1.A[end], v1_leg.A[end]; rtol=rtol)
    @test isapprox(v1.Q[end], v1_leg.Q[end]; rtol=rtol)
    @test isapprox(v1.u[end], v1_leg.u[end]; rtol=rtol)
    @test isapprox(v2.A[end], v2_leg.A[end]; rtol=rtol)
    @test isapprox(v2.Q[end], v2_leg.Q[end]; rtol=rtol)
    @test isapprox(v2.u[end], v2_leg.u[end]; rtol=rtol)
    @test isapprox(v3.A[1],   v3_leg.A[1];   rtol=rtol)
    @test isapprox(v3.Q[1],   v3_leg.Q[1];   rtol=rtol)
    @test isapprox(v3.u[1],   v3_leg.u[1];   rtol=rtol)
end
