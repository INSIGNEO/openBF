#=
Copyright 2020 INSIGNEO Institute for in silico Medicine

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


"""
    loadSimulationFiles(input_filename :: String)

Load, check, and return `.yml` input file content.
"""
function loadSimulationFiles(input_filename :: String)
    data = loadYAMLFile(input_filename)

    checkInputFile(data)

    return data
end


"""
    loadYAMLFile(filename :: String)

Check existence and open a `.yml` file and return the content as `Dict{Any,Any}`.
"""
function loadYAMLFile(filename :: String)
    if ~isfile(filename)
        error("missing file $filename")
    end

    return YAML.load(open(filename))
end


"""
    checkInputFile(input_filename :: String)

Check YAML input file. Raise errors if any of the important variables is ill
defined.
"""
function checkInputFile(data :: Dict{Any,Any})
    checkSections(data)
    checkNetwork(data["network"])
end


"""
    checkSections(data :: Dict{Any,Any})

Look for the four sections in the input data. Run integrity checks for `blood` and
`solver` sections.
"""
function checkSections(data :: Dict{Any,Any})
    keys = ["project name", "network", "blood", "solver"]
    for key in keys
        if ~haskey(data, key)
            error("missing section $key in YAML input file")
        end
    end

    checkSection(data, "blood", ["mu", "rho"])
    checkSection(data, "solver", ["Ccfl", "cycles", "convergence tolerance"])

    if ~haskey(data["solver"], "jump")
        data["solver"]["jump"] = 100
    end
end


"""
    checkSection(data :: Dict{Any,Any}, section :: String, keys :: Array{String,1})

Look for a list of keys in the given data section.
"""
function checkSection(data :: Dict{Any,Any}, section :: String, keys :: Array{String,1})
    for key in keys
        if ~haskey(data[section], key)
            error("missing $key in $section section")
        end
    end
end


"""
    checkNetwork(network :: Array{Dict{Any,Any},1})

Loop trough the network and run checks on each single vessel. Check also if at least one
inlet and one oulet has been defined.
"""
function checkNetwork(network :: Array{Dict{Any,Any},1})
    has_inlet = false
    has_outlet = false
    nodes = Dict{Int,Int}()
    for i = 1:length(network)
        checkVessel(i, network[i])

        if haskey(network[i], "inlet")
            has_inlet = true
        end
        if haskey(network[i], "outlet")
            has_outlet = true
        end

        # check max number of vessels per node
        if ~haskey(nodes, network[i]["sn"])
            nodes[network[i]["sn"]] = 1
        else
            nodes[network[i]["sn"]] += 1
        end
        if ~haskey(nodes, network[i]["tn"])
            nodes[network[i]["tn"]] = 1
        else
            nodes[network[i]["tn"]] += 1
        end
        if nodes[network[i]["sn"]] > 3
            error("too many vessels connected at node $(network[i]["sn"])")
        elseif nodes[network[i]["tn"]] > 3
            error("too many vessels connected at node $(network[i]["tn"])")
        end
    end

    # outlet nodes must be defined and all network connnected
    for i = 1:length(network)
        if nodes[network[i]["tn"]] == 1
            if ~haskey(network[i], "outlet")
                error("outlet not defined for vessel $(network[i]["label"]),
                    check connectivity")
            end
        end
        delete!(nodes, network[i]["tn"])
    end
    for k in keys(network)
        if k == 1
            continue
        else
            if network[k] == 1
                error("node $(k) is not an outlet and it is not connected
                 to at least two vessels, check connectivity")
            end
        end
    end

    if ~has_inlet
        error("missing inlet(s) definition")
    end

    if ~has_outlet
        error("missing outlet(s) definition")
    end
end


"""
    checkVessel(i :: Int, vessel :: Dict{Any,Any})

Check if all the important parameters are defined for the current vessel. If the vessel
has been indicated to be an inlet/outlet segment, check also for the definition of
boundary condition parameters.
"""
function checkVessel(i :: Int, vessel :: Dict{Any,Any})
    keys = ["label", "sn", "tn", "L", "E"]
    for key in keys
        if ~haskey(vessel, key)
            error("vessel $i is missing $key value")
        end
    end

    if vessel["sn"] == vessel["tn"]
        error("vessel $i has same sn and tn")
    end

    if ~haskey(vessel, "R0")
        if ~haskey(vessel, "Rp") && ~haskey(vessel, "Rd")
            error("vessel $i is missing lumen radius value(s)")
        end
    else
        if vessel["R0"] > 0.05
            @warn "$(vessel["label"]) radius larger than 5cm!"
        end
    end

    if haskey(vessel, "inlet")
        if ~haskey(vessel, "inlet file")
            error("inlet vessel $i is missing the inlet file path")
        elseif ~isfile(vessel["inlet file"])
            file_path = vessel["inlet file"]
            error("vessel $i inlet file $file_path not found")
        end

        if ~haskey(vessel, "inlet number")
            error("inlet vessel $i is missing the inlet number")
        end
    end

    if haskey(vessel, "outlet")
        outlet = vessel["outlet"]
        if outlet == "wk3"
            if ~haskey(vessel, "R1") || ~haskey(vessel, "Cc")
                error("outlet vessel $i is missing three-element windkessel values")
            end
        elseif outlet == "wk2"
            if ~haskey(vessel, "R1") || ~haskey(vessel, "Cc")
                error("outlet vessel $i is missing two-element windkessel values")
            end
        elseif outlet == "reflection"
            if ~haskey(vessel, "Rt")
                error("outlet vessel $i is missing reflection coefficient value")
            end
        end
    end
end


"""
    makeResultsFolder(data :: Dict{Any,Any})

Create results folder and cd in.
"""
function makeResultsFolder(data :: Dict{Any,Any}, input_filename :: String)
    project_name = data["project name"]

    if ~haskey(data, "results folder")
        r_folder = join([project_name, "_results"])
    else
        r_folder = data["results folder"]
    end

    if isdir(r_folder) == false
      mkdir(r_folder)
    end

    cp(input_filename, r_folder*"/"*input_filename)

    cd(r_folder)
end


"""
    buildBlood(project_constants :: Dict{Any,Any})

Return an instance of type `::Blood`.
"""
function buildBlood(blood_data :: Dict{Any,Any})
    mu = blood_data["mu"]
    rho = blood_data["rho"]
    rho_inv = 1.0/rho

    return Blood(mu, rho, rho_inv)
end


"""
    buildArterialNetwork(network :: Array{Dict{Any,Any},1}, blood :: Blood, jump :: Int)

Populate `vessels` and `edges` lists containing mechanical properties and topology og the
system, respectively.
"""
function buildArterialNetwork(network :: Array{Dict{Any,Any},1}, blood :: Blood,
                              jump :: Int)
    vessels = [buildVessel(1, network[1], blood, jump)]
    edges = zeros(Int, length(network), 3)
    @inbounds edges[1,1] = vessels[1].ID
    @inbounds edges[1,2] = vessels[1].sn
    @inbounds edges[1,3] = vessels[1].tn

    for i = 2:length(network)
        push!(vessels, buildVessel(i, network[i], blood, jump))
        @inbounds edges[i,1] = vessels[i].ID
        @inbounds edges[i,2] = vessels[i].sn
        @inbounds edges[i,3] = vessels[i].tn
    end

    return vessels, edges
end


"""
    buildVessel(vessel_data :: Dict{Any,Any}, blood :: Blood, jump :: Int)

Allocate arrays, set initial conditions, and create a new instance of `Vessel` with data
from the `.yml`.
"""
function buildVessel(ID :: Int, vessel_data :: Dict{Any,Any}, blood :: Blood, jump :: Int)
    vessel_name = vessel_data["label"]
    sn = vessel_data["sn"]
    tn = vessel_data["tn"]
    L = vessel_data["L"]
    E = vessel_data["E"]

    Rp, Rd = computeRadii(vessel_data)
    Pext = getPext(vessel_data)
    M, dx, invDx, halfDx, invDxSq = meshVessel(vessel_data, L)
    h0 = initialiseThickness(vessel_data, M)
    outlet, Rt, R1, R2, Cc = addOutlet(vessel_data)
    viscT = computeViscousTerm(vessel_data, blood)
    inlet, heart = buildHeart(vessel_data)
    phi = getPhi(vessel_data)

    Q = zeros(Float64, M)
    P = zeros(Float64, M)
    A = zeros(Float64, M)
    u = zeros(Float64, M)
    c = zeros(Float64, M)
    A0 = zeros(Float64, M)
    R0 = zeros(Float64, M)
    s_A0 = zeros(Float64, M)
    beta = zeros(Float64, M)
    vA = zeros(Float64, M+2)
    vQ = zeros(Float64, M+2)
    Al = zeros(Float64, M+2)
    Ar = zeros(Float64, M+2)
    Ql = zeros(Float64, M+2)
    Qr = zeros(Float64, M+2)
    wallE = zeros(Float64, M)
    gamma = zeros(Float64, M)
    slope = zeros(Float64, M)
    wallVa = zeros(Float64, M)
    wallVb = zeros(Float64, M)
    inv_A0 = zeros(Float64, M)
    dU = zeros(Float64, 2, M+2)
    Fl = zeros(Float64, 2, M+2)
    Fr = zeros(Float64, 2, M+2)
    s_inv_A0 = zeros(Float64, M)
    slopesA = zeros(Float64, M+2)
    slopesQ = zeros(Float64, M+2)
    Q_t = zeros(Float64, jump, 6)
    P_t = zeros(Float64, jump, 6)
    A_t = zeros(Float64, jump, 6)
    u_t = zeros(Float64, jump, 6)
    c_t = zeros(Float64, jump, 6)
    Q_l = zeros(Float64, jump, 6)
    P_l = zeros(Float64, jump, 6)
    A_l = zeros(Float64, jump, 6)
    u_l = zeros(Float64, jump, 6)
    c_l = zeros(Float64, jump, 6)
    flux  = zeros(Float64, 2, M+2)
    uStar = zeros(Float64, 2, M+2)
    s_15_gamma = zeros(Float64, M)
    gamma_ghost = zeros(Float64, M+2)

    s_pi = sqrt(pi)
    s_pi_E_over_sigma_squared = s_pi*E/0.75
    one_over_rho_s_p = 1.0/(3.0*blood.rho*s_pi)
    radius_slope = computeRadiusSlope(Rp, Rd, L)

    ah = 0.2802
    bh = -5.053e2
    ch = 0.1324
    dh = -0.1114e2

    if h0 == 0.0
        Rmean = 0.5*(Rp + Rd)
        h0 = computeThickness(Rmean, ah, bh, ch, dh)
    end
    Cv = 0.5*s_pi*phi*h0/(blood.rho*0.75)

    @fastmath @inbounds for i = 1:M
      R0[i] = radius_slope*(i - 1)*dx + Rp
      A0[i] = pi*R0[i]*R0[i]
      s_A0[i] = sqrt(A0[i])
      inv_A0[i] = 1.0/A0[i]
      s_inv_A0[i] = sqrt(inv_A0[i])
      A[i] = A0[i]
      beta[i] = s_inv_A0[i]*h0*s_pi_E_over_sigma_squared
      gamma[i] = beta[i]*one_over_rho_s_p/R0[i]
      s_15_gamma[i] = sqrt(1.5*gamma[i])
      gamma_ghost[i+1] = gamma[i]
      P[i] = pressure(1.0, beta[i], Pext)
      c[i] = waveSpeed(A[i], gamma[i])
      wallE[i] = 3.0*beta[i]*radius_slope*inv_A0[i]*s_pi*blood.rho_inv
      if phi != 0.0
          wallVb[i] = Cv*s_inv_A0[i]*invDxSq
          wallVa[i] = 0.5*wallVb[i]
      end
    end

    gamma_ghost[1] = gamma[1]
    gamma_ghost[end] = gamma[end]

    if outlet == "wk2"
        R1, R2 = computeWindkesselInletImpedance(R2, blood, A0, gamma)
        outlet = "wk3"
    end

    U00A = A0[1]
    U01A = A0[2]
    UM1A = A0[M]
    UM2A = A0[M-1]

    U00Q = 0.0
    U01Q = 0.0
    UM1Q = 0.0
    UM2Q = 0.0

    W1M0 = u[end] - 4.0*c[end]
    W2M0 = u[end] + 4.0*c[end]

    node2 = convert(Int, floor(M*0.25))
    node3 = convert(Int, floor(M*0.5))
    node4 = convert(Int, floor(M*0.75))

    Pcn = 0.0

    last_A_name = join((vessel_name,"_A.last"))
    last_Q_name = join((vessel_name,"_Q.last"))
    last_u_name = join((vessel_name,"_u.last"))
    last_c_name = join((vessel_name,"_c.last"))
    last_P_name = join((vessel_name,"_P.last"))

    out_A_name = join((vessel_name,"_A.out"))
    out_Q_name = join((vessel_name,"_Q.out"))
    out_u_name = join((vessel_name,"_u.out"))
    out_c_name = join((vessel_name,"_c.out"))
    out_P_name = join((vessel_name,"_P.out"))

    return Vessel(vessel_name, ID, sn, tn, inlet, heart,
                  M, dx, invDx, halfDx,
                  beta, gamma, s_15_gamma, gamma_ghost,
                  A0, s_A0, inv_A0, s_inv_A0, Pext,
                  viscT, wallE, wallVa, wallVb,
                  A, Q, u, c, P,
                  A_t, Q_t, u_t, c_t, P_t,
                  A_l, Q_l, u_l, c_l, P_l,
                  W1M0, W2M0,
                  U00A, U00Q, U01A, U01Q, UM1A, UM1Q, UM2A, UM2Q,
                  last_P_name, last_Q_name, last_A_name,
                  last_c_name, last_u_name,
                  out_P_name, out_Q_name, out_A_name,
                  out_c_name, out_u_name,
                  node2, node3, node4,
                  Rt, R1, R2, Cc,
                  Pcn,
                  slope, flux, uStar, vA, vQ,
                  dU, slopesA, slopesQ,
                  Al, Ar, Ql, Qr, Fl, Fr,
                  outlet)
end


"""
    computeRadiusSlope(Rd :: Float64, Rp :: Float64, L :: Float64).

Calculate the slope for the lumen radius linear tapering as

(Rd - Rp)/L

"""
function computeRadiusSlope(Rp :: Float64, Rd :: Float64, L :: Float64)
    return (Rd - Rp)/L
end


"""
    computeThickness(R0i :: Float64)

Compute wall thickness based on local lumen radius as

h_{0i} = R_{0i} ( a_h e^{b_h R_{0i}} + c_h e^{d_h R_{0i}} )

where

- ah = 0.2802
- bh = -5.053e2 m^{-1}
- ch = 0.1324
- dh = -0.1114e2 m^{-1}

as reported in

> Avolio AP. Multi-branched model of the human arterial system. Medical and Biological Engineering and Computing. 1980 Nov 1;18(6):709-18.
"""
function computeThickness(R0i :: Float64,
                          ah :: Float64, bh :: Float64, ch :: Float64, dh :: Float64)
    return R0i*(ah*exp(bh*R0i) + ch*exp(dh*R0i))
end


"""
    computeRadii(vessel_data :: Dict{Any,Any})

If only a constant lumen radius is defined, return the same value for proximal and
distal variables, `Rp` and `Rd`, respectively.
"""
function computeRadii(vessel :: Dict{Any,Any})
    if ~haskey(vessel, "R0")
        Rp = vessel["Rp"]
        Rd = vessel["Rd"]
        return Rp, Rd
    else
        R0 = vessel["R0"]
        return R0, R0
    end
end


"""
    getPext(vessel :: Dict{Any,Any})

Extract Pext value for current vessels; return default `Pext = 0.0` if no value is
specified.
"""
function getPext(vessel :: Dict{Any,Any})
    if ~haskey(vessel, "Pext")
        return 0.0
    else
        return vessel["Pext"]
    end
end


"""
    getPhi(vessel :: Dict{Any,Any})

Extract `phi` value for current vessels; return default `phi = 0.0` if no value is
specified.
"""
function getPhi(vessel :: Dict{Any,Any})
    if ~haskey(vessel, "phi")
        return 0.0
    else
        return vessel["phi"]
    end
end


"""
    meshVessel(vessel :: Dict{Any,Any}, L :: Float64)

Pre-compute Delta x, frac{1}{Delta x} and frac{Delta x}{2} for the current
vessel. The Delta x is computed as

Delta x = frac{L}{M}

where `M` is the maximum value between `5` (the minimum needed by the solver), the value
defined in the `.yml`, and `ceil(L*1e3)` (which would make Delta x = 1mm).
"""
function meshVessel(vessel :: Dict{Any,Any}, L :: Float64)
    if ~haskey(vessel, "M")
        M = maximum([5, convert(Int, ceil(L*1e3))])
    else
        M = vessel["M"]
        M = maximum([5, M, convert(Int, ceil(L*1e3))])
    end

    dx = L/M
    invDx = M/L
    halfDx = 0.5*dx
    invDxSq = invDx*invDx

    return M, dx, invDx, halfDx, invDxSq
end


"""
    initialiseThickness(vessel :: Dict{Any,Any}, M :: Int)

If vessel thickness is not specified in the `.yml` file, return a `zeroes` array to be
filed with radius dependent values. Otherwise, return an array with the specified `h0`.
"""
function initialiseThickness(vessel :: Dict{Any,Any}, M :: Int)
    if ~haskey(vessel, "h0")
        return 0.0
    else
        return vessel["h0"]
    end
end


"""
    addOutlet(vessel :: Dict{Any,Any})

Parse outlet information for the current vessel and return windkessel and reflection
coeffiecient values.
"""
function addOutlet(vessel :: Dict{Any,Any})
    if haskey(vessel, "outlet")
        outlet = vessel["outlet"]
        if outlet == "wk3"
            Rt = 0.0
            R1 = vessel["R1"]
            R2 = vessel["R2"]
            Cc = vessel["Cc"]
        elseif outlet == "wk2"
            Rt = 0.0
            R1 = 0.0
            R2 = vessel["R1"]
            Cc = vessel["Cc"]
        elseif outlet == "reflection"
            Rt = vessel["Rt"]
            R1 = 0.0
            R2 = 0.0
            Cc = 0.0
        end
    else
        outlet = "none"
        Rt = 0.0
        R1 = 0.0
        R2 = 0.0
        Cc = 0.0
    end

    return outlet, Rt, R1, R2, Cc
end


"""
    computeViscousTerm(vessel_data :: Dict{Any,Any}, blood :: Blood)

Return

2(gamma_v + 2)pi*mu

where gamma_v (`gamma_profile`) is either specified in the vessel definition or
assumed equal to `9` (plug-flow).
"""
function computeViscousTerm(vessel_data :: Dict{Any,Any}, blood :: Blood)
    if haskey(vessel_data, "gamma_profile")
        gamma_profile = vessel_data["gamma_profile"]
    else
        gamma_profile = 9
    end
    return 2*(gamma_profile + 2)*pi*blood.mu*blood.rho_inv
end


"""
    buildHeart(vessel_data)

If the current vessel is an inlet vessel, return `true` flag and an `Heart` struct.
"""
function buildHeart(vessel :: Dict{Any,Any})
    if haskey(vessel, "inlet")
        inlet_type = vessel["inlet"]
        input_data = loadInletData(vessel["inlet file"])
        cardiac_period = input_data[end, 1]
        inlet_number = vessel["inlet number"]

        return true, Heart(inlet_type, cardiac_period, input_data, inlet_number)
    else
        return false, Heart("none", 0.0, zeros(1,2), 0)
    end
end


"""
    loadInletData(inlet_file :: String)

Read discretised inlet data from inlet file. Return an `Array{Float64, 2}` whose first
columun contains time variable and second column contains the inlet time function.
"""
function loadInletData(inlet_file :: String)
    return readdlm(inlet_file)
end


"""
    computeWindkesselInletImpedance(R2 :: Float64, blood :: Blood, A0 :: Array{Float64,1},
                                    gamma :: Array{Float64,1})

In case only one peripheral resistance is defined (two-element windkessel), the second
one is set as equal to the outlet vessel impedance.
"""
function computeWindkesselInletImpedance(R2 :: Float64, blood :: Blood,
    A0 :: Array{Float64,1}, gamma :: Array{Float64,1})

    R1 = blood.rho*waveSpeed(A0[end], gamma[end])/A0[end]
    R2 -= R1

    return R1, R2
end


# http://carlobaldassi.github.io/ArgParse.jl/stable/index.html
function parseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "input_filename"
            help = ".yml input file name"
            required = true
        "--verbose", "-v"
            help = "Print STDOUT - default false"
            action = :store_true
        "--out_files", "-f"
            help = "Save complete results story rather than only the last cardiac cycle"
            action = :store_true
    end

    return parse_args(s)
end
