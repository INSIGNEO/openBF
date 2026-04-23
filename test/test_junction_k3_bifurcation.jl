# Validate solve_junction! at k=3 against the existing join_vessels!(v1,v2,v3)
# bifurcation solver. Bifurcation uses static pressure (use_total_pressure=false).

@testset "Junction solve k=3 (bifurcation)" begin
    blood = openBF.Blood(Dict("rho" => 1060.0, "mu" => 0.004))
    cfg_parent = Dict{Any,Any}(
        "label" => "parent", "sn" => 1, "tn" => 2,
        "L" => 0.086, "E" => 500_000.0, "R0" => 7.581e-3, "h0" => 0.9e-3,
        "gamma_profile" => 9,
    )
    cfg_d1 = Dict{Any,Any}(
        "label" => "d1", "sn" => 2, "tn" => 3,
        "L" => 0.085, "E" => 700_000.0, "R0" => 5.492e-3, "h0" => 0.68e-3,
        "gamma_profile" => 9,
        "R1" => 6.8123e7, "R2" => 3.1013e9, "Cc" => 3.6664e-10,
        "inlet_impedance_matching" => false,
    )
    cfg_d2 = Dict{Any,Any}(
        "label" => "d2", "sn" => 2, "tn" => 4,
        "L" => 0.085, "E" => 700_000.0, "R0" => 5.492e-3, "h0" => 0.68e-3,
        "gamma_profile" => 9,
        "R1" => 6.8123e7, "R2" => 3.1013e9, "Cc" => 3.6664e-10,
        "inlet_impedance_matching" => false,
    )

    v1 = openBF.Vessel(cfg_parent, blood, 100, ["P", "Q", "A", "u"])
    v2 = openBF.Vessel(cfg_d1,    blood, 100, ["P", "Q", "A", "u"])
    v3 = openBF.Vessel(cfg_d2,    blood, 100, ["P", "Q", "A", "u"])

    # Perturb boundary faces to a non-trivial state.
    v1.A[end] *= 1.03
    v1.Q[end]  = v1.A[end] * 0.30
    v1.u[end]  = v1.Q[end] / v1.A[end]

    v2.A[1]   *= 0.97
    v2.Q[1]    = v2.A[1] * 0.14
    v2.u[1]    = v2.Q[1] / v2.A[1]

    v3.A[1]   *= 0.98
    v3.Q[1]    = v3.A[1] * 0.15
    v3.u[1]    = v3.Q[1] / v3.A[1]

    v1_leg = deepcopy(v1)
    v2_leg = deepcopy(v2)
    v3_leg = deepcopy(v3)

    # Legacy bifurcation solve.
    openBF.join_vessels!(v1_leg, v2_leg, v3_leg)

    # Generic junction solve (static pressure).
    jc = Junction(2, [1, 2, 3], [:outlet, :inlet, :inlet])
    solve_junction!(jc, [v1, v2, v3], blood)

    rtol = 1e-10

    @test isapprox(v1.A[end], v1_leg.A[end]; rtol=rtol)
    @test isapprox(v1.Q[end], v1_leg.Q[end]; rtol=rtol)
    @test isapprox(v1.u[end], v1_leg.u[end]; rtol=rtol)
    @test isapprox(v2.A[1],   v2_leg.A[1];   rtol=rtol)
    @test isapprox(v2.Q[1],   v2_leg.Q[1];   rtol=rtol)
    @test isapprox(v2.u[1],   v2_leg.u[1];   rtol=rtol)
    @test isapprox(v3.A[1],   v3_leg.A[1];   rtol=rtol)
    @test isapprox(v3.Q[1],   v3_leg.Q[1];   rtol=rtol)
    @test isapprox(v3.u[1],   v3_leg.u[1];   rtol=rtol)
end
