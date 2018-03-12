#=
Copyright 2017 INSIGNEO Institute for in silico Medicine

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

# http://carlobaldassi.github.io/ArgParse.jl/stable/index.html
function parseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "project_name"
            help = "Project name"
            required = true
        "--verbose", "-v"
            help = "Print STDOUT - default false"
            action = :store_true
        "--clean", "-c"
            help = "Clean input, .out, and .temp files at the end of the simulation - default false"
            action = :store_true
    end

    return parse_args(s)
end


# check basic files
function checkInputFiles(project_name :: String)
    f_const = join([project_name, "_constants.yml"])
    f_inlet = join([project_name, "_inlet.dat"])
    f_model = join([project_name, ".csv"])

    # check if all the input files are in the current folder
    for f in [f_const, f_inlet, f_model]
        if isfile(f) == false
            error("$f file missing")
        end
    end
end

# check additional inlets
function checkInletFiles(project_name :: String, number_of_inlets :: Int64,
                         inlets :: Array{String, 1})

    for inlet_idx = 2:number_of_inlets
        f_additional_inlet = join([project_name, "_", inlet_idx, "_inlet.dat"])

        if isfile(f_additional_inlet) == false
            error("number_of_inlets = $number_of_inlets, but $f file is missing")
        end

        push!(inlets, f_additional_inlet)
    end

    return inlets
end


# create new folder for results and copy input files
function copyInputFilesToResultsFolder(project_name :: String, inlets :: Array{String, 1})
    # make results folder
    r_folder = join([project_name, "_results"])
    if isdir(r_folder) == false
      mkdir(r_folder)
    end

    # copy input files in results folder
    f_const = join([project_name, "_constants.yml"])
    cp(f_const, join([r_folder, "/", f_const]), remove_destination=true)

    f_model = join([project_name, ".csv"])
    cp(f_model, join([r_folder, "/", f_model]), remove_destination=true)

    for f_inlet in inlets
        cp(f_inlet, join([r_folder, "/", f_inlet]), remove_destination=true)
    end

    cd(r_folder)
end


function loadInletData(inlet_file :: String)
    return readdlm(inlet_file)
end


function buildHeart(project_constants :: Dict{Any,Any}, inlet_data :: Array{Float64,2},
                    inlet_number :: Int64)

    cardiac_period = inlet_data[end, 1]
    initial_flow = 0.0

    return Heart(project_constants["inlet_type"], cardiac_period, inlet_data,
                 initial_flow, inlet_number)
end


function buildHearts(project_constants :: Dict{Any,Any}, inlet_names :: Array{String, 1})

    hearts = []
    for i = 1:length(inlet_names)
        inlet_file = inlet_names[i]
        inlet_data = loadInletData(inlet_file)
        push!(hearts, buildHeart(project_constants, inlet_data, i))
    end

    return hearts
end


function buildBlood(project_constants :: Dict{Any,Any})
    mu = project_constants["mu"]
    rho = project_constants["rho"]
    gamma_profile = project_constants["gamma_profile"]

    nu = mu/rho
    Cf = 8*pi*nu
    rho_inv = 1/rho
    viscT = 2*(gamma_profile + 2)*pi*mu

    return Blood(mu, rho, rho_inv, Cf, gamma_profile, viscT)
end


function checkConstants(project_constants :: Dict{Any,Any})

    fundamental_parameters = ["inlet_type", "mu", "rho", "number_of_inlets"]
    for key in fundamental_parameters
        if ~haskey(project_constants, key)
            error("$key not defined in <project>.yml")
        end
    end

    not_so_important_parameters = ["gamma_profile", "Ccfl", "cycles", "initial_pressure"]
    default_values = [9, 0.9, 100, 0.0]
    for i = 1:length(not_so_important_parameters)
        key = not_so_important_parameters[i]
        if ~haskey(project_constants, key)
            default_value = default_values[i]
            warn("$key not defined in <project>.yml, assuming $default_value")
            project_constants[key] = default_values[i]
        end
    end

    return project_constants
end


function loadConstants(project_constants_file :: String)
    return YAML.load(open(project_constants_file))
end


function loadSimulationFiles(project_name :: String)

    checkInputFiles(project_name)

    f_inlet = join([project_name, "_inlet.dat"])
    inlets = [f_inlet]

    # load constants
    f_const = join([project_name, "_constants.yml"])
    project_constants = loadConstants(f_const)
    project_constants = checkConstants(project_constants)

    # check for additional inlets
    number_of_inlets = project_constants["number_of_inlets"]
    if number_of_inlets > 1
        inlets = checkInletFiles(project_name, number_of_inlets, inlets)
    end

    # load inlets data
    hearts = buildHearts(project_constants, inlets)

    # load blood data
    blood = buildBlood(project_constants)

    # estimate total simulation time
    total_time = project_constants["cycles"]*hearts[1].cardiac_T

    # make results folder and copy input files
    copyInputFilesToResultsFolder(project_name, inlets)

    f_model = join([project_name, ".csv"])
    model, model_header = readModelData(f_model)

    return [project_constants, model, hearts, blood, total_time]
end


function parseModelRow(model_row :: Array{Any,1})
    vessel_name = model_row[1]

    if typeof(vessel_name) != SubString{String}
        error("vessel $vessel_name: the vessel name must be a string beginning with a literal (a,..., z, A,..., Z)")
    else
        vessel_name = convert(String, vessel_name)
    end

    sn = convert(Int64, model_row[2])
    tn = convert(Int64, model_row[3])
    rn = convert(Int64, model_row[4])
    L = convert(Float64, model_row[5])
    M = convert(Int64, model_row[6])
    Rp = convert(Float64, model_row[7])
    Rd = convert(Float64, model_row[8])
    E = convert(Float64, model_row[9])
    Pext = convert(Float64, model_row[10])

    if model_row[11] == ""
        Rt = ""
        R1 = ""
        R2 = ""
        Cc = ""
    elseif model_row[11] != "" && model_row[12] == ""
        Rt = convert(Float64, model_row[11])
        R1 = ""
        R2 = ""
        Cc = ""
    elseif model_row[11] != "" && model_row[12] != "" && model_row[13] == ""
        Rt = ""
        R1 = ""
        R2 = convert(Float64, model_row[11])
        Cc = convert(Float64, model_row[12])
    elseif model_row[11] != "" && model_row[12] != "" && model_row[13] != ""
        Rt = ""
        R1 = convert(Float64, model_row[11])
        R2 = convert(Float64, model_row[12])
        Cc = convert(Float64, model_row[13])
    end

    return vessel_name, sn, tn, rn, L, M, Rp, Rd, E, Pext, [Rt, R1, R2, Cc]
end

function meshVessel(L :: Float64, M :: Int64)
    M = maximum([5, M, convert(Int, ceil(L*1e3))])
    dx = L/M
    invDx = M/L
    halfDx = 0.5*dx

    return dx, invDx, halfDx
end

# no outlet BC
function checkCapillaries(BCout :: Array{String,1}, blood :: Blood,
                          A0 :: Array{Float64,1}, gamma :: Array{Float64,1})
    BCout = [0.0, 0.0, 0.0, 0.0]
    return BCout, "none"
end

# parse outlet BCs
function checkCapillaries(BCout :: Array{Any,1}, blood :: Blood,
                          A0 :: Array{Float64,1}, gamma :: Array{Float64,1})
    if BCout[1] != "" && BCout[2] == ""
        BCout[2:4] = [0.0, 0.0, 0.0]
        return BCout, "reflection"

    elseif BCout[1] == "" && BCout[2] == "" && BCout[3] != ""
        BCout[1] = 0.0
        BCout[2] = blood.rho*waveSpeed(A0[end], gamma[end])/A0[end]
        BCout[3] -= BCout[2]
        return BCout, "wk3"

    elseif BCout[1] == "" && BCout[2] != ""
        BCout[1] = 0.0
        return BCout, "wk3"
    end
end


function buildArterialNetwork(model :: Array{Any, 2}, heart :: Heart, blood :: Blood)
    vessels = [buildVessel(1, model[1,:], heart, blood)]
    edge_list = zeros(Int64, size(model)[1], 4)
    edge_list[1,1] = vessels[1].ID
    edge_list[1,2] = vessels[1].sn
    edge_list[1,3] = vessels[1].tn
    edge_list[1,4] = vessels[1].inlet_idx

    for i = 2:size(model)[1]
        push!(vessels, buildVessel(i, model[i,:], heart, blood))
        edge_list[i,1] = vessels[i].ID
        edge_list[i,2] = vessels[i].sn
        edge_list[i,3] = vessels[i].tn
        edge_list[i,4] = vessels[i].inlet_idx
    end

    return vessels, edge_list
end


function buildVessel(ID :: Int64, model_row :: Array{Any,1},
                     heart :: Heart, blood :: Blood)

    vessel_name, sn, tn, rn, L, M, Rp, Rd, E, Pext, BCout = parseModelRow(model_row)

    dx, invDx, halfDx = meshVessel(L, M)

    Q = zeros(Float64, M)
    P = zeros(Float64, M)
    A = zeros(Float64, M)
    u = zeros(Float64, M)
    c = zeros(Float64, M)
    A0 = zeros(Float64, M)
    R0 = zeros(Float64, M)
    h0 = zeros(Float64, M)
    beta = zeros(Float64, M)
    vA = zeros(Float64, M+2)
    vQ = zeros(Float64, M+2)
    Al = zeros(Float64, M+2)
    Ar = zeros(Float64, M+2)
    Ql = zeros(Float64, M+2)
    Qr = zeros(Float64, M+2)
    gamma = zeros(Float64, M)
    dA0dx = zeros(Float64, M)
    slope = zeros(Float64, M)
    dTaudx = zeros(Float64, M)
    inv_A0 = zeros(Float64, M)
    dU = zeros(Float64, 2, M+2)
    Fl = zeros(Float64, 2, M+2)
    Fr = zeros(Float64, 2, M+2)
    s_inv_A0 = zeros(Float64, M)
    slopesA = zeros(Float64, M+2)
    slopesQ = zeros(Float64, M+2)
    flux  = zeros(Float64, 2, M+2)
    uStar = zeros(Float64, 2, M+2)
    gamma_ghost = zeros(Float64, M+2)
    half_beta_dA0dx = zeros(Float64, M)

    s_pi = sqrt(pi)
    s_pi_E_over_sigma_squared = s_pi*E/0.75
    one_over_rho_s_p = 1/(3*blood.rho*s_pi)
    radius_slope = (Rd-Rp)/(M-1)
    ah = 0.2802
    bh = -5.053e2
    ch = 0.1324
    dh = -0.1114e2

    for i = 1:M
      R0[i] = radius_slope*(i - 1)*dx + Rp
      h0[i] = R0[i]*(ah*exp(bh*R0[i]) + ch*exp(dh*R0[i]))
      A0[i] = pi*R0[i]*R0[i]
      A[i] = A0[i]
      Q[i] = heart.initial_flow
      u[i] = Q[i]/A[i]
      inv_A0[i] = 1/A0[i]
      s_inv_A0[i] = sqrt(inv_A0[i])
      dA0dx[i] = 2*pi*R0[i]*radius_slope
      dTaudx[i] = sqrt(pi)*E*radius_slope*1.3*(h0[i]/R0[i] + R0[i]*(ah*bh*exp(bh*R0[i]) + ch*dh*exp(dh*R0[i])))
      beta[i] = s_inv_A0[i]*h0[i]*s_pi_E_over_sigma_squared
      gamma[i] = beta[i]*one_over_rho_s_p/R0[i]
      gamma_ghost[i+1] = gamma[i]
      half_beta_dA0dx[i] = beta[i]*0.5*dA0dx[i]
      P[i] = pressure(A[i], A0[i], beta[i], Pext)
      c[i] = waveSpeed(A[i], gamma[i])
    end

    gamma_ghost[1] = gamma[1]
    gamma_ghost[end] = gamma[end]

    BCout, outlet = checkCapillaries(BCout, blood, A0, gamma)
    Rt = BCout[1]
    R1 = BCout[2]
    R2 = BCout[3]
    Cc = BCout[4]

    U00A = A0[1]
    U01A = A0[2]
    UM1A = A0[M]
    UM2A = A0[M-1]

    U00Q = heart.initial_flow
    U01Q = heart.initial_flow
    UM1Q = heart.initial_flow
    UM2Q = heart.initial_flow

    W1M0 = u[end] - 4*c[end]
    W2M0 = u[end] + 4*c[end]

    node2 = convert(Int, floor(M*0.25))
    node3 = convert(Int, floor(M*0.5))
    node4 = convert(Int, floor(M*0.75))

    Pcn = 0.

    temp_A_name = join((vessel_name,"_A.temp"))
    temp_Q_name = join((vessel_name,"_Q.temp"))
    temp_u_name = join((vessel_name,"_u.temp"))
    temp_c_name = join((vessel_name,"_c.temp"))
    temp_P_name = join((vessel_name,"_P.temp"))

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

    temp_A = open(temp_A_name, "w")
    temp_Q = open(temp_Q_name, "w")
    temp_u = open(temp_u_name, "w")
    temp_c = open(temp_c_name, "w")
    temp_P = open(temp_P_name, "w")

    last_A = open(last_A_name, "w")
    last_Q = open(last_Q_name, "w")
    last_u = open(last_u_name, "w")
    last_c = open(last_c_name, "w")
    last_P = open(last_P_name, "w")

    out_A = open(out_A_name, "w")
    out_Q = open(out_Q_name, "w")
    out_u = open(out_u_name, "w")
    out_c = open(out_c_name, "w")
    out_P = open(out_P_name, "w")

    close(last_A)
    close(last_Q)
    close(last_u)
    close(last_c)
    close(last_P)

    close(out_A)
    close(out_Q)
    close(out_u)
    close(out_c)
    close(out_P)

    return Vessel(vessel_name, ID, sn, tn, rn, M,
                  dx, invDx, halfDx,
                  beta, gamma, gamma_ghost, half_beta_dA0dx,
                  A0, inv_A0, s_inv_A0, dA0dx, dTaudx, Pext,
                  A, Q, u, c, P,
                  W1M0, W2M0,
                  U00A, U00Q, U01A, U01Q, UM1A, UM1Q, UM2A, UM2Q,
                  temp_P_name, temp_Q_name, temp_A_name,
                  temp_c_name, temp_u_name,
                  last_P_name, last_Q_name, last_A_name,
                  last_c_name, last_u_name,
                  out_P_name, out_Q_name, out_A_name,
                  out_c_name, out_u_name,
                  temp_P, temp_Q, temp_A, temp_c, temp_u,
                  last_P, last_Q, last_A, last_c, last_u,
                  node2, node3, node4,
                  Rt, R1, R2, Cc,
                  Pcn,
                  slope, flux, uStar, vA, vQ,
                  dU, slopesA, slopesQ,
                  Al, Ar, Ql, Qr, Fl, Fr,
                  outlet)
end

# All the pre-processing and initialisation functions are herein contained.



# *function* __`readModelData`__ $\rightarrow$ `model_matrix::Array{Any, 2}`
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `model_data`     `.csv` file containing model data (see
#                  [tutorial.jl](../index.html#project)).
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning:
# ----------------------------------------------------------------------------
# `model_data.csv` is parsed by `readdlm` by using `,` as a separator. The
# first header row is discarded (`[2:end]`) as well as the last column
# (`[1:end-1]`) because it contains
# only blank spaces. The remaining is stored in matrix `m` and returned.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ---------------- -----------------------------------------------------------
# `m`              `::Array{Any, 2}` model matrix.
# ----------------------------------------------------------------------------
# <a name="readModelData"></a>

function readModelData(model_csv :: String)

  m, h = readdlm(model_csv, ',', header=true)

  for i = 1:length(h)
      h[i] = strip(h[i])
  end

  return m, h
end
