# Validate solve_junction! at k=2 against the existing join_vessels! conjunction solver.
# Conjunction uses total pressure (static + dynamic), so use_total_pressure=true.

@testset "Junction solve k=2 (conjunction)" begin
    blood = openBF.Blood(Dict("rho" => 1060.0, "mu" => 0.004))
    cfg_parent = Dict{Any,Any}(
        "label" => "parent", "sn" => 1, "tn" => 2,
        "L" => 0.1, "E" => 500_000.0, "R0" => 0.003, "h0" => 3e-4,
        "gamma_profile" => 9,
    )
    cfg_child = Dict{Any,Any}(
        "label" => "child", "sn" => 2, "tn" => 3,
        "L" => 0.1, "E" => 700_000.0, "R0" => 0.002, "h0" => 2e-4,
        "gamma_profile" => 9,
        "R1" => 1e8, "R2" => 1e9, "Cc" => 1e-10,
        "inlet_impedance_matching" => false,
    )

    # Construct vessels, perturb boundary faces to a non-trivial state.
    v1 = openBF.Vessel(cfg_parent, blood, 100, ["P", "Q", "A", "u"])
    v2 = openBF.Vessel(cfg_child,  blood, 100, ["P", "Q", "A", "u"])

    v1.A[end] *= 1.02
    v1.Q[end]  = v1.A[end] * 0.25
    v1.u[end]  = v1.Q[end] / v1.A[end]

    v2.A[1]   *= 0.98
    v2.Q[1]    = v2.A[1] * 0.24
    v2.u[1]    = v2.Q[1] / v2.A[1]

    # Deep-copy for the legacy solver so both start from identical state.
    v1_leg = deepcopy(v1)
    v2_leg = deepcopy(v2)

    # Legacy conjunction solve.
    openBF.join_vessels!(v1_leg, v2_leg, blood.rho)

    # Generic junction solve (total pressure, matching conjunction).
    jc = Junction(2, [1, 2], [:outlet, :inlet]; use_total_pressure=true)
    solve_junction!(jc, [v1, v2], blood)

    rtol = 1e-10

    @test isapprox(v1.A[end], v1_leg.A[end]; rtol=rtol)
    @test isapprox(v1.Q[end], v1_leg.Q[end]; rtol=rtol)
    @test isapprox(v1.u[end], v1_leg.u[end]; rtol=rtol)
    @test isapprox(v2.A[1],   v2_leg.A[1];   rtol=rtol)
    @test isapprox(v2.Q[1],   v2_leg.Q[1];   rtol=rtol)
    @test isapprox(v2.u[1],   v2_leg.u[1];   rtol=rtol)
end
