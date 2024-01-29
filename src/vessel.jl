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

struct Blood
    mu::Float64
    rho::Float64
    Cf::Float64
end

function Blood(config::Dict)
    μ = config["mu"]
    ρ = config["rho"]
    Cf = 8π * μ / ρ
    Blood(μ, ρ, Cf)
end

mutable struct Vessel
    label::SubString{String}

    #Topological notation
    sn::Int64
    tn::Int64

    #Numerical constants
    M::Int64
    dx::Float64
    invDx::Float64
    halfDx::Float64
    Ccfl::Float64

    #Physical constants
    beta::Vector{Float64}
    Cv::Vector{Float64}
    gamma::Vector{Float64}
    gamma_ghost::Vector{Float64}
    A0::Vector{Float64}
    dA0dx::Vector{Float64}
    dTaudx::Vector{Float64}
    Pext::Float64
    gamma_profile::Int64

    #Iterative solution
    A::Vector{Float64}
    Q::Vector{Float64}
    u::Vector{Float64}
    c::Vector{Float64}
    P::Vector{Float64}

    #Riemann invariants
    W1M0::Float64
    W2M0::Float64

    #Ghost cells
    U00A::Float64
    U00Q::Float64
    U01A::Float64
    U01Q::Float64

    UM1A::Float64
    UM1Q::Float64
    UM2A::Float64
    UM2Q::Float64

    #Saving locations
    node2::Int64
    node3::Int64
    node4::Int64

    #Peripheral boundary condition
    Rt::Float64

    R1::Float64
    R2::Float64
    Cc::Float64
    Pc::Float64

    #Slope
    slope::Vector{Float64}

    #MUSCLArrays
    flux::Array{Float64,2}
    uStar::Array{Float64,2}

    vA::Vector{Float64}
    vQ::Vector{Float64}

    dU::Array{Float64,2}

    slopesA::Vector{Float64}
    slopesQ::Vector{Float64}

    Al::Vector{Float64}
    Ar::Vector{Float64}

    Ql::Vector{Float64}
    Qr::Vector{Float64}

    Fl::Array{Float64,2}
    Fr::Array{Float64,2}

end

wave_speed(A::Float64, gamma::Float64) = sqrt(1.5 * gamma * sqrt(A))

pressure(A::Float64, A0::Float64, beta::Float64, Pext::Float64) =
    Pext + beta * (sqrt(A / A0) - 1.0)

function mesh(config::Dict{Any,Any})
    L = config["L"]
    M = get(config, "M", 5)
    M = max(M, 5)
    M = max(M, ceil(Int, config["L"] * 1e-3))

    dx = L / M
    invDx = M / L
    halfDx = 0.5 * dx

    L, M, dx, invDx, halfDx
end

function radii(config::Dict{Any,Any})
    ~haskey(config, "R0") &&
        ~haskey(config, "Rp") &&
        error("missing radius in $(config["label"])")
    R0 = get(config, "R0", 0.0)
    Rp = get(config, "Rp", R0)
    Rd = get(config, "Rd", Rp)
    Rp, Rd
end

function Vessel(config::Dict{Any,Any}, b::Blood, Ccfl::Float64)
    vessel_name = config["label"]
    sn = config["sn"]
    tn = config["tn"]
    Rp, Rd = radii(config)
    L, M, dx, invDx, halfDx = mesh(config)
    E = config["E"]

    Pext = get(config, "Pext", 0.0)
    initial_pressure = get(config, "initial_pressure", 0.0)
    initial_flow = get(config, "initial_flow", 0.0)

    sigma = 0.5

    A0 = zeros(Float64, M)
    R0 = zeros(Float64, M)
    h0 = zeros(Float64, M) .+ get(config, "h0", 0.0)

    dA0dx = zeros(Float64, M)
    dTaudx = zeros(Float64, M)
    radius_slope = (Rd - Rp) / (M - 1)
    ah = 0.2802
    bh = -5.053e2
    ch = 0.1324
    dh = -0.1114e2
    for i = 1:M
        R0[i] = radius_slope * (i - 1) * dx + Rp
        if ~haskey(config, "h0")
            h0[i] = R0[i] * (ah * exp(bh * R0[i]) + ch * exp(dh * R0[i]))
        end
        A0[i] = pi * R0[i] * R0[i]
        dA0dx[i] = 2 * pi * R0[i] * radius_slope
        dTaudx[i] =
            sqrt(pi) *
            E *
            radius_slope *
            1.3 *
            (
                h0[i] / R0[i] +
                R0[i] * (ah * bh * exp(bh * R0[i]) + ch * dh * exp(dh * R0[i]))
            )
    end

    beta = sqrt.(pi ./ A0) .* h0 * E / (1 - sigma^2)
    # beta = zeros(Float64, M) .+ get(config, "beta", 0.0)

    gamma = beta ./ (3 * b.rho * R0 * sqrt(pi))

    phi_v = get(config, "phi_v", 0.0) # default elastic
    Cv = 2*sqrt(pi)*phi_v.*h0./(b.rho.*sqrt.(A0))

    gamma_ghost = zeros(Float64, M + 2)
    gamma_ghost[2:M+1] = gamma
    gamma_ghost[1] = gamma[1]
    gamma_ghost[end] = gamma[end]

    A = zeros(Float64, M) + A0
    Q = zeros(Float64, M) .+ initial_flow
    u = zeros(Float64, M) + Q ./ A
    c = wave_speed.(A, gamma)
    P = pressure.(A, A0, beta, Pext)

    U00A = A0[1]
    U01A = A0[2]
    UM1A = A0[M]
    UM2A = A0[M-1]

    U00Q = initial_flow
    U01Q = initial_flow
    UM1Q = initial_flow
    UM2Q = initial_flow

    W1M0 = u[end] - 4 * c[end]
    W2M0 = u[end] + 4 * c[end]

    node2 = round(Int, M * 0.25)
    node3 = round(Int, M * 0.5)
    node4 = round(Int, M * 0.75)

    Rt = get(config, "Rt", 0.0)
    R1 = get(config, "R1", 0.0)
    R2 = get(config, "R2", 0.0)
    Cc = get(config, "Cc", 0.0)
    if R2 == 0.0
        R2 = R1 - b.rho * wave_speed(A0[end], gamma[end]) / A0[end]
    end

    # `Pc` is the pressure through the peripheral compliance of the three
    # elements windkessel. It is set to zero to simulate the pressure at the
    # artery-vein interface.
    Pc = 0.0

    #Slope
    slope = zeros(Float64, M)

    # MUSCL arrays
    flux = zeros(Float64, 2, M + 2)
    uStar = zeros(Float64, 2, M + 2)

    vA = zeros(Float64, M + 2)
    vQ = zeros(Float64, M + 2)

    dU = zeros(Float64, 2, M + 2)

    slopesA = zeros(Float64, M + 2)
    slopesQ = zeros(Float64, M + 2)

    Al = zeros(Float64, M + 2)
    Ar = zeros(Float64, M + 2)

    Ql = zeros(Float64, M + 2)
    Qr = zeros(Float64, M + 2)

    Fl = zeros(Float64, 2, M + 2)
    Fr = zeros(Float64, 2, M + 2)

    gamma_profile = get(config, "gamma_profile", 2)

    Vessel(
        vessel_name,
        sn,
        tn,
        M,
        dx,
        invDx,
        halfDx,
        Ccfl,
        beta,
        Cv,
        gamma,
        gamma_ghost,
        A0,
        dA0dx,
        dTaudx,
        Pext,
        gamma_profile,
        A,
        Q,
        u,
        c,
        P,
        W1M0,
        W2M0,
        U00A,
        U00Q,
        U01A,
        U01Q,
        UM1A,
        UM1Q,
        UM2A,
        UM2Q,
        node2,
        node3,
        node4,
        Rt,
        R1,
        R2,
        Cc,
        Pc,
        slope,
        flux,
        uStar,
        vA,
        vQ,
        dU,
        slopesA,
        slopesQ,
        Al,
        Ar,
        Ql,
        Qr,
        Fl,
        Fr,
    )
end
