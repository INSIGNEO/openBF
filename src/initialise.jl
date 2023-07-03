#=
Copyright 2022 INSIGNEO Institute for in silico Medicine

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
    data
end


"""
    loadYAMLFile(filename :: String)

Check existence and open a `.yml` file and return the content as `Dict{Any,Any}`.
"""
function loadYAMLFile(filename :: String)
    ~isfile(filename) && error("missing file $filename")
    YAML.load(open(filename))
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
raiseMissingKey(data::Dict{Any,Any}, key::String) = ~haskey(data, key) && error("missing section $key in YAML input file")

function checkSections(data :: Dict{Any,Any})
    keys = ["project name", "network", "blood", "solver"]
    raiseMissingKey.(Ref(data), keys)
    checkSection.(Ref(data), "blood", ["mu", "rho"])
    checkSection.(Ref(data), "solver", ["Ccfl", "cycles", "convergence tolerance"])
    ~haskey(data["solver"], "jump") && (data["solver"]["jump"] = 100)
end


"""
    checkSection(data :: Dict{Any,Any}, section :: String, keys :: Array{String,1})

Look for a list of keys in the given data section.
"""
checkSection(d::Dict{Any,Any}, s::String, k::String) = ~haskey(d[s], k) && error("missing $k in $s")


"""
    checkNetwork(network :: Array{Dict{Any,Any},1})

Loop trough the network and run checks on each single vessel. Check also if at least one
inlet and one oulet has been defined.
"""
function checkNetwork(network :: Array{Dict{Any,Any},1})
    has_inlet = false
    has_outlet = false
    nodes = Dict{Int,Int}()
    for i = eachindex(network)
        checkVessel(i, network[i])

        haskey(network[i], "inlet") && (has_inlet = true)
        haskey(network[i], "outlet") && (has_outlet = true)

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

    # outlet nodes must be defined
    for i = 1:length(network)
        nodes[network[i]["tn"]] == 1 && ~haskey(network[i], "outlet") && error("outlet not @ $(network[i]["label"])")
    end

    ~has_inlet && error("missing inlet(s) definition")
    ~has_outlet && error("missing outlet(s) definition")
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

    vessel["sn"] == vessel["tn"] && error("vessel $i has same sn and tn")

    if ~haskey(vessel, "R0")
        ~haskey(vessel, "Rp") && ~haskey(vessel, "Rd") && error("vessel $i is missing lumen radius value(s)")
    else
        vessel["R0"] > 0.05 && @warn "$(vessel["label"]) radius larger than 5cm!"
    end

    if haskey(vessel, "inlet")
        if ~haskey(vessel, "inlet file")
            error("inlet vessel $i is missing the inlet file path")
        elseif ~isfile(vessel["inlet file"])
            file_path = vessel["inlet file"]
            error("vessel $i inlet file $file_path not found")
        end

        ~haskey(vessel, "inlet number") && error("inlet vessel $i is missing the inlet number")
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

    # delete existing folder and results!
    if isdir(r_folder)
        rm(r_folder, recursive=true)
    end
    mkdir(r_folder)

    cp(input_filename, r_folder*"/"*input_filename, force=true)
    copyInletFiles(data, r_folder)

    cd(r_folder)
end


"""
    copyInletFiles(data :: Dict{Any,Any}, r_folder :: String)

Copy inlet .dat files to results folder.
"""
function copyInletFiles(data :: Dict{Any,Any}, r_folder :: String)
    for vessel in data["network"]
        if haskey(vessel, "inlet file")
            cp(vessel["inlet file"], r_folder*"/"*vessel["inlet file"], force=true)
        end
    end
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
    n = length(network)+1
    grafo = Graphs.SimpleDiGraph(n)
    v = buildVessel(1, network[1], blood, jump)
    vessels = Dict(Graphs.Edge(v.sn, v.tn)=>v)
    Graphs.add_edge!(grafo, v.sn, v.tn)

    for i = 2:length(network)
        v = buildVessel(i, network[i], blood, jump)
        Graphs.add_edge!(grafo, v.sn, v.tn)
        vessels[Graphs.Edge(v.sn, v.tn)] = v
    end

    edges = Graphs.edges(grafo)

    # vessels = [buildVessel(1, network[1], blood, jump)]
    # edges = zeros(Int, length(network), 3)
    # edges[1,1] = vessels[1].ID
    # edges[1,2] = vessels[1].sn
    # edges[1,3] = vessels[1].tn

    # for i = 2:length(network)
    #     push!(vessels, buildVessel(i, network[i], blood, jump))
    #     edges[i,1] = vessels[i].ID
    #     edges[i,2] = vessels[i].sn
    #     edges[i,3] = vessels[i].tn
    # end

    return grafo, vessels, edges
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
    outlet, Rt, R1, R2, Cc = addOutlet(vessel_data)
    viscT = computeViscousTerm(vessel_data, blood) # ---------------> ???
    inlet, heart = buildHeart(vessel_data)

    ah = 0.2802
    bh = -5.053e2
    ch = 0.1324
    dh = -0.1114e2
    h0 = zeros(Float64, M)
    R0 = zeros(Float64, M)
    radius_slope = computeRadiusSlope(Rp, Rd, L)
    dTaudx = zeros(Float64, M)

    for i = 1:M
        R0[i] = radius_slope*(i - 1)*dx + Rp
        h0[i] = computeThickness(R0[i], ah, bh, ch, dh)
        dTaudx[i] = sqrt(pi)*E*radius_slope*1.3*(h0[i]/R0[i]+R0[i]*(ah*bh*exp(bh*R0[i]) + ch*dh*exp(dh*R0[i])))
    end
    A0 = pi.*R0.^2
    s_A0 = sqrt.(A0)
    inv_A0 = 1.0 ./ A0
    s_inv_A0 = sqrt.(inv_A0)
    dA0dx = 2*pi.*R0.*radius_slope
    
    sigma = 0.5
    beta  = sqrt.(pi./A0) .* h0*E / (1 - sigma^2)

    gamma = beta ./ (3*blood.rho*R0*sqrt(pi))
    s_15_gamma = sqrt.(1.5 .*gamma)

    gamma_ghost = zeros(Float64, M+2)
    gamma_ghost[2:M+1] = gamma
    gamma_ghost[1] = gamma[1]
    gamma_ghost[end] = gamma[end]

    A = zeros(Float64, M)  + A0
    P = pressure.( A, A0, beta, Pext)
    Q = zeros(Float64, M)
    u = zeros(Float64, M)  + Q./A
    c = waveSpeed.(A, gamma)

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

    vA = zeros(Float64, M+2)
    vQ = zeros(Float64, M+2)
    Al = zeros(Float64, M+2)
    Ar = zeros(Float64, M+2)
    Ql = zeros(Float64, M+2)
    Qr = zeros(Float64, M+2)
    wallE = zeros(Float64, M)  # ???
    slope = zeros(Float64, M)
    wallVa = zeros(Float64, M) # ???
    wallVb = zeros(Float64, M) # ???
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
    flux = zeros(Float64, 2, M+2)
    uStar = zeros(Float64, 2, M+2)
    dU = zeros(Float64, 2, M+2)
    Fl = zeros(Float64, 2, M+2)
    Fr = zeros(Float64, 2, M+2)
    s_inv_A0 = zeros(Float64, M)
    slopesA = zeros(Float64, M+2)
    slopesQ = zeros(Float64, M+2)
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
                  outlet, dTaudx, dA0dx)
end


"""
    computeRadiusSlope(Rd :: Float64, Rp :: Float64, L :: Float64).

Calculate the slope for the lumen radius linear tapering as

(Rd - Rp)/L

"""
computeRadiusSlope(Rp::Float64, Rd::Float64, M::Float64) = (Rd - Rp)/(M-1)


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
computeThickness(R0i::Float64, ah::Float64, bh::Float64, ch::Float64, dh::Float64) = (exp(bh*R0i)ah + exp(dh*R0i)ch)R0i


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
    end
    R0 = vessel["R0"]
    return R0, R0
end


"""
    getPext(vessel :: Dict{Any,Any})

Extract Pext value for current vessels; return default `Pext = 0.0` if no value is
specified.
"""
getPext(vessel :: Dict{Any,Any}) = ~haskey(vessel, "Pext") ? 10.0e3 : vessel["Pext"]


"""
    getPhi(vessel :: Dict{Any,Any})

Extract `phi` value for current vessels; return default `phi = 0.0` if no value is
specified.
"""
getPhi(vessel :: Dict{Any,Any}) = ~haskey(vessel, "phi") ? 0.0 : vessel["phi"]


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
        M = maximum([5, M]) #, convert(Int, ceil(L*1e3))])
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
initialiseThickness(vessel::Dict{Any,Any}, M::Int) = ~haskey(vessel, "h0") ? 0.0 : vessel["h0"]


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
    gamma_profile = haskey(vessel_data, "gamma_profile") ? vessel_data["gamma_profile"] : 9
    2(gamma_profile + 2)pi*blood.mu #*blood.rho_inv
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
    end
    return false, Heart("none", 0.0, zeros(1,2), 0)
end


"""
    loadInletData(inlet_file :: String)

Read discretised inlet data from inlet file. Return an `Array{Float64, 2}` whose first
columun contains time variable and second column contains the inlet time function.
"""
loadInletData(inlet_file :: String) = readdlm(inlet_file)


"""
    computeWindkesselInletImpedance(R2 :: Float64, blood :: Blood, A0 :: Vector{Float64},
                                    gamma :: Vector{Float64})

In case only one peripheral resistance is defined (two-element windkessel), the second
one is set as equal to the outlet vessel impedance.
"""
function computeWindkesselInletImpedance(R2 :: Float64, blood :: Blood,
    A0 :: Vector{Float64}, gamma :: Vector{Float64})

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
            help = "Print STDOUT (default false)"
            action = :store_true
        "--out_files", "-f"
            help = "Save complete results story rather than only the last cardiac cycle (default true)"
            action = :store_true
        "--conv_ceil", "-c"
            help = "Ceil convergence value to 100 mmHg (default true)"
            action = :store_true
    end

    return parse_args(s)
end
