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

# All the pre-processing and initialisation functions are herein contained.

# *function* __`projectPreamble`__
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `project_name`   `::String` project name.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning:
# ----------------------------------------------------------------------------
# This function checks whether all the project files are present in the
# working directory. Then, the directory structure to save the output files
# is created. Eventually, `openBF` logo and version information are printed to
# screen.
# ----------------------------------------------------------------------------
# <a name="projectPreamble"></a>
function projectPreamble(project_name :: String, no_out :: Bool, no_inputs :: Bool, number_of_inlets :: Int64)

  # `project_constants.jl`, `project_inlet.dat`, and `project.csv` file names
  # strings are created with `project_name` variable.
  p_const = join([project_name, "_constants.jl"])
  p_inlet = join([project_name, "_inlet.dat"])

  if number_of_inlets > 1
      p_inlets = []
      for inlet_idx = 2:number_of_inlets
          push!(p_inlets, join([project_name, "_", inlet_idx, "_inlet.dat"]))
      end
  end

  p_model = join([project_name, ".csv"])
  v_model = join([project_name, "_veins.csv"])
  # The existence of these three files is checked in the working directory.
  # If a file is missing an error is raised.
  if (isfile(p_const)) == false
    error("$p_const file missing")
  end

  if (isfile(p_inlet)) == false
    println("$p_inlet file missing")
  end

  if (isfile(p_model)) == false
    error("$p_model file missing")
  end
  # `r_folder` is the string used to create the directory containing the
  # results.
  r_folder = join([project_name, "_results"])
  # If the results directory does not exist, it is created.
  if (isdir(r_folder)) == false
    mkdir(r_folder)
  end
  # Project files are moved to the results directory.
  cp(p_const, join([r_folder, "/", p_const]), remove_destination=true)
  cp(p_inlet, join([r_folder, "/", p_inlet]), remove_destination=true)

  if number_of_inlets > 1
     for i = 1:number_of_inlets-1
         cp(p_inlets[i], join([r_folder, "/", p_inlets[i]]), remove_destination=true)
     end
  end

  cp(p_model, join([r_folder, "/", p_model]), remove_destination=true)
  if (isfile(v_model)) == true
    cp(v_model, join([r_folder, "/", v_model]), remove_destination=true)
  end
  # All files in the working directory are deleted except those containing
  # the `project_name` string by the
  # [`cleanLibrary`](initialise.html#cleanLibrary) function.
  # cleanLibrary(project_name)
  # run(`rm clean_lib.sh`)
  # The working directory is changed to the results directory.
  cd(r_folder)
  # Bash scripts to handle I/O are written and saved in the working directory
  # by [`writeScripts`](initialise.html#writeScripts) function.
  # writeScripts(no_out, no_inputs)

end

# Two bash scripts are written in two `.sh` files.
#
# *   `appender.sh` is used to concatenate two files which names should be
#     given as script arguments. It can be used as
#
#         sh appender.sh file1 file2
#
# *   `cleaner.sh` is used at the end of the simulation to remove all the
#     temporary files from the working directory.
#
# These `.sh` are saved in the working directory and used through the `julia`
# command `run`.

# *function* __`writeScripts`__
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# Two bash scripts are written in two `.sh` files.
# ----------------------------------------------------------------------------
#
# *   `appender.sh` is used to concatenate two files which names should be
#     given as script arguments. It can be used as
#
#         sh appender.sh file1 file2
#
# *   `cleaner.sh` is used at the end of the simulation to remove all the
#     temporary files from the working directory.
#
# These `.sh` are saved in the working directory and used through the `julia`
# command `run`.
# <a name="writeScripts"></a>
function writeScripts(no_out :: Bool, no_inputs :: Bool)

  appsh = open("appender.sh", "w")
  write(appsh, "#!/bin/bash", "\n")
  write(appsh, "cat \$1 >> \$2")
  close(appsh)

  # clesh = open("cleaner.sh", "w")
  # write(clesh, "#!/bin/bash", "\n")
  # write(clesh, "rm *.temp\n")
  # write(clesh, "rm appender.sh")
  # if no_out == true
  #     write(clesh, "\nrm *.out")
  # end
  # if no_inputs == true
  #     write(clesh, "\nrm *.csv")
  #     write(clesh, "\nrm *.dat")
  #     write(clesh, "\nrm *.jl")
  # end
  # close(clesh)

end

# *function* __`cleanLibrary`__
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `project_name`   `::String` project name.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning:
# ----------------------------------------------------------------------------
# This function writes and executes a bash scripts which delete all the
# files with the name not containing the string `project_name`.
# ----------------------------------------------------------------------------
# <a name="cleanLibrary"></a>
function cleanLibrary(project_name)

  clibsh = open("clean_lib.sh", "w")
  write(clibsh, "#!/bin/bash", "\n")
  write(clibsh, "find . -type f -not -name '*$project_name*' | xargs rm -rf", "\n")
  write(clibsh, "rm -rf doc/")
  close(clibsh)

  #run(`sh clean_lib.sh`)

end

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
function readModelData(model_data)

  m = readdlm(model_data, ',')
  m = m[2:end, 1:end-1]

  return m
end

# *function* __`initialiseVessel`__ $\rightarrow$
# `vessel::`[`Vessel`]( html#Vessel)
#
# ----------------------------------------------------------------------------
# Parameters:
# ------------------- --------------------------------------------------------
# `m`                 `::Array{Any, 2}` model matrix row.
#
# `ID`                `::Int` topological ID of the vessel in the arterial
#                     tree.
#
# `h`                 [`::Heart`]( html#Heart) data structure
#                     containing inlet boundary conditions.
#
# `b`                 [`::Blood`]( html#Blood) data structure
#                     containing blood mechanical properties.
#
# `Pext`              `::Float` external pressure.
#
# `initial_pressure`  `::Float` initial pressure to be set along the entire
#                     vessel.
#
# `Ccfl`              `::Float` Courant-Friedrichs-Lewy number to [compute
#                     local $\Delta t$](godunov.html#calculateDeltaT).
# ----------------------------------------------------------------------------
# <a name="initialiseVessel"></a>
function initialiseVessel(m :: Array{Any, 1}, ID :: Int64, h :: Heart,
                              b :: Blood,
                              initial_pressure :: Float64,
                              Ccfl :: Float64)

  # --------------------------------------------------------------------------
  # Functioning:
  # --------------------------------------------------------------------------
  # Model matrix row `m` is parsed  to retrieve vessel geometrical and
  # mechanical properties; conversion to `::Int` type is operated where
  # needed.
  # --------------------------------------------------------------------------
  vessel_name = m[1]
  sn = convert(Int, m[2])
  tn = convert(Int, m[3])
  rn = convert(Int, m[4])
  L = m[5]

  #M = convert(Int, m[6])
  M = maximum([5, convert(Int, m[6]), convert(Int, ceil(L*1e3))])

  Rp = m[7]
  Rd = m[8]

#   h0          =     m[7]

  E = m[9]
  Pext = convert(Float64, m[10])
  # Poisson's ratio `sigma` is set by default to 0.5 because vessel wall is
  # assumed to be incompressible.
  sigma = 0.5
  # $\Delta x$ for local numerical discretisation: lenght/number of nodes.
  # <a name="dx"></a>
  dx = L/M
  invDx = M/L
  halfDx = 0.5*dx

  # Unstressed cross-sectional area `A0` is computed from `R0`.
#   A0 = pi*R0*R0

#   Rp = 4e-3
#   Rd = 2e-3

#   Rp = r0
#   Rd = r0

  A0 = zeros(Float64, M)
  inv_A0 = zeros(Float64, M)
  s_inv_A0 = zeros(Float64, M)
  R0 = zeros(Float64, M)
  h0 = zeros(Float64, M)

  # Elastic constants for trans-mural [pressure](converter.html#pressure)
  # (`beta`) and [wave speed](converter.html#waveSpeed) (`gamma`)
  # calculation:
  # $$
  #   \beta = \sqrt{\frac{\pi}{A_0}} \frac{h_0 E}{1 - \sigma^2}, \quad
  #   \gamma = \frac{\beta}{3 \rho R_0 \sqrt{\pi}} .
  # $$
  beta = zeros(Float64, M)
  gamma = zeros(Float64, M)

  dA0dx = zeros(Float64, M)
  dTaudx = zeros(Float64, M)
  radius_slope = (Rd-Rp)/(M-1)
  ah = 0.2802
  bh = -5.053e2
  ch = 0.1324
  dh = -0.1114e2

  half_beta_dA0dx = zeros(Float64, M)

  for i = 1:M
    R0[i] = radius_slope*(i-1)*dx + Rp
    h0[i] = R0[i] * (ah * exp(bh*R0[i]) + ch * exp(dh*R0[i])) #1e-3
    A0[i] = pi*R0[i]*R0[i]
    inv_A0[i] = 1/A0[i]
    s_inv_A0[i] = sqrt(inv_A0[i])
    dA0dx[i] = 2*pi*R0[i]*radius_slope
    dTaudx[i] = sqrt(pi)*E*radius_slope*1.3*(h0[i]/R0[i]+R0[i]*(ah*bh*exp(bh*R0[i]) + ch*dh*exp(dh*R0[i])))

    beta[i] = sqrt(pi/A0[i])*h0[i]*E/(1 - sigma*sigma)
    gamma[i] = beta[i]/(3*b.rho*R0[i]*sqrt(pi))

    half_beta_dA0dx[i] = beta[i]*0.5*dA0dx[i]
  end

  gamma_ghost = zeros(Float64, M+2)
  gamma_ghost[2:M+1] = gamma
  gamma_ghost[1] = gamma[1]
  gamma_ghost[end] = gamma[end]

  # Initialise conservative (`A` and `Q`) and primitive (`u`, `c`, and `P`)
  # variables arrays with initial values. `A` is taken as the unstressed
  # cross-sectional area `A0` and `Q` is everywhere the `initial_flow`
  # specified in [`project_constants.jl`](../index.html#project_constants).
  # The longitudinal velocity `u` comes from the definition of `Q`
  # $$
  #     Q = u A .
  # $$
  # `c` and `P` are computed with `beta` and `gamma` by
  # [`pressure`](converter.html#pressure) and
  # [`waveSpeed`](converter.html#waveSpeed) functions.
  A = zeros(Float64, M)  + A0
  Q = zeros(Float64, M)  + h.initial_flow
  u = zeros(Float64, M)  + Q./A
  c = zeros(Float64, M)
  c = waveSpeed(A, gamma, c)
  P = zeros(Float64, M)
  P = pressure( A, A0, beta, Pext, P)
  # [Ghost cells](godunov.html) are initialised
  # as the internal nodes.
  U00A = A0[1]
  U01A = A0[2]
  UM1A = A0[M]
  UM2A = A0[M-1]

  U00Q = h.initial_flow
  U01Q = h.initial_flow
  UM1Q = h.initial_flow
  UM2Q = h.initial_flow
  # Forward and backward
  # [Riemann invariants](converter.html#riemannInvariants)
  # at the vessel outlet node.
  W1M0 = u[end] - 4*c[end]
  W2M0 = u[end] + 4*c[end]
  # Temporary files names are built by `join`ing the vessel `label` and a
  # string referring to the quantity stored, and a `.temp` extension. Same
  # procedure is followed for results files with `.out` extension. Data input
  # and output are handled by [`IOutils.jl`](IOutils.html) functions.
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
  # Temporary data and output files are created with `open` command. Data from
  # temporary files to output files are
  # [transferred](IOutils.jl#transferTempToOut) by means of a bash script.
  # Thus, the `IOStream` to `.out` files must be `close`d. Writing in `.temp`
  # is handled within julia, hence `.temp` are left open, ready to be
  # written.
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

  node2 = convert(Int, floor(M*0.25))
  node3 = convert(Int, floor(M*0.5))
  node4 = convert(Int, floor(M*0.75))

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
  # The outlet boundary condition is specified by the last columns in the
  # `project.csv` file.
  #
  # 9 columns means that only the reflection coefficient `Rt` has been
  # specified. Since the three elements windkessel parameters must be
  # specified anyway, they are set to zero.
  if length(m) == 11

    Rt = m[11]

    R1 = 0.
    R2 = 0.
    Cc = 0.
  # 10 columns means that the three element windkessel is specified with
  # Reymond (2009) notation. Peripheral resistance (`R2`) is supplied along
  # with the outlet impedance (`R1`). To decouple these two parameters, `R1`
  # is computed as the vessel characteristic impedance
  # $$
  #   R_1 = Z_c = \rho \frac{c}{A_0}.
  # $$
  elseif length(m) == 12

    R1 = b.rho*waveSpeed(A0[end], gamma[end])/A0[end]
    #R1 = 1e5 #b.rho * 13.3/(A0*(2*R0*1e3)^0.3)
    R2 = m[11] - R1
    Cc = m[12]

    Rt = 0.
  # 11 columns means that the three parameters are specified.
  elseif length(m) >= 13

    R1 = m[11]
    R2 = m[12]
    Cc = m[13]

    Rt = 0.

  end

#     elseif length(m) == 10

#     R1 = b.rho*waveSpeed(A0[end], gamma[end])/A0[end]
#     #R1 = 1e5 #b.rho * 13.3/(A0*(2*R0*1e3)^0.3)
#     R2 = m[ 9] - R1
#     Cc = m[10]

#     Rt = 0.
#   # 11 columns means that the three parameters are specified.
#   elseif length(m) >= 11

#     R1 = b.rho*waveSpeed(A0[end], gamma[end])/A0[end] #m[ 9]
#     R2 = m[ 9] +m[10]#- R1 #m[10]
#     Cc = m[11]

#     Rt = 0.

#   end


  # `Pcn` is the pressure through the peripheral compliance of the three
  # elements windkessel. It is set to zero to simulate the pressure at the
  # artery-vein interface.
  Pcn = 0.

  #Slope
  slope = zeros(Float64, M)

  # MUSCL arrays
  flux  = zeros(Float64, 2, M+2)
  uStar = zeros(Float64, 2, M+2)

  vA = zeros(Float64, M+2)
  vQ = zeros(Float64, M+2)

  dU = zeros(Float64, 2, M+2)

  slopesA = zeros(Float64, M+2)
  slopesQ = zeros(Float64, M+2)

  Al = zeros(Float64, M+2)
  Ar = zeros(Float64, M+2)

  Ql = zeros(Float64, M+2)
  Qr = zeros(Float64, M+2)

  Fl = zeros(Float64, 2, M+2)
  Fr = zeros(Float64, 2, M+2)
#
# ----------------------------------------------------------------------------
# Returns:
# -------------- -------------------------------------------------------------
# `vessel_data`  [`::Vessel`]( html#Vessel) data structure
#                containing vessel mechanical and geometrical properties.
# ----------------------------------------------------------------------------
  vessel_data =  Vessel(vessel_name,
                    ID, sn, tn, rn,
                    M,
                    dx, invDx, halfDx,
                    Ccfl,
                    beta, gamma, gamma_ghost, half_beta_dA0dx,
#                     R0,
                    A0, inv_A0, s_inv_A0, dA0dx, dTaudx, Pext,
                    A, Q,
                    u, c, P,
                    W1M0, W2M0,
                    U00A, U00Q, U01A, U01Q,
                    UM1A, UM1Q, UM2A, UM2Q,
                    temp_P_name, temp_Q_name, temp_A_name,
                    temp_c_name, temp_u_name,
                    last_P_name, last_Q_name, last_A_name,
                    last_c_name, last_u_name,
                    out_P_name, out_Q_name, out_A_name,
                    out_c_name, out_u_name,
                    temp_P, temp_Q, temp_A, temp_c, temp_u,
                    last_P, last_Q, last_A, last_c, last_u,
                    node2, node3, node4,
                    Rt,
                    R1, R2, Cc,
                    Pcn,
                    slope,
                    flux, uStar, vA, vQ,
                    dU, slopesA, slopesQ,
                    Al, Ar, Ql, Qr, Fl, Fr)

  return vessel_data
end

# *function* __`loadInputData`__ $\rightarrow$ `::Array{Float, 2}`
#
# ----------------------------------------------------------------------------
# Parameters:
# ------------------- --------------------------------------------------------
# `project_name`      `project` name string.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# `project_inlet.dat` is parsed by `readdlm` function and its content is
# saved in a `in_data` matrix. Its first column contains time and the second
# column contains flow waveform data.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# --------- ------------------------------------------------------------------
# `in_data` `::Array{Float, 2}` `[time, flow]` matrix for inlet boundary
# condition.
# ----------------------------------------------------------------------------
# <a name="loadInputData"></a>
function loadInputData(project_name)
  in_data = readdlm(join([project_name, "_inlet.dat"]))
  return in_data
end

# *function* __`loadGlobalConstants`__ $\rightarrow$ `::Heart`,
# `::Blood`, `::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# --------------------- ------------------------------------------------------
# `project_name`        `project` name string.
#
# `inlet_BC_switch`     `::Int` inlet boundary condtion selector:
#
#                        1 $\rightarrow$ heaviside sine function
#
#                        2 $\rightarrow$ gaussian
#
#                        3 $\rightarrow$ from `project_inlet.dat`
#
# `cycles`              `::Int` number of cardiac cycles to be simulated.
#
# `rho`                 `::Float` blood density
#
# `mu`                  `::Float` blood dynamic viscosity
# ----------------------------------------------------------------------------
#
# --------------------------------------------------------------------------
# Functioning
# --------------------------------------------------------------------------
# This function load into [heart]( html#Heart) and
# [blood]( html#Blood) datastructures the parameters constant
# within all vessels.
# --------------------------------------------------------------------------
# <a name="loadGlobalConstants"></a>
function loadGlobalConstants(project_name,
                             inlet_BC_switch :: Int64,
                             inlet_type :: String,
                             cycles :: Int64,
                             rho :: Float64, mu:: Float64, gamma_profile :: Int64,
                             number_of_inlets :: Int64)
  # `heart` data structure is filled depending on the inlet boundary condition
  # chosen by the user. When the inlet flow function is given
  # (`inlet_BC_swithc`=3), the flow waveform is imported by
  # [`loadInputData`](initialise.html#loadInputData) function. The
  # `cardiac_period` is taken as the last row in the time column of the input
  # file. `sys_T`, `initial_flow`, and `flow_amplitude` quantities will not
  # be used and they are set to zero.
  if inlet_BC_switch == 3
    input_data = loadInputData(project_name)

#     if input_data[1,2] < 1e-6
#       error("Inlet flow initial value < 1e-6")
#     end

    cardiac_period = input_data[end, 1]
    sys_T = 0.
    initial_flow = 0.
    flow_amplitude = 0.
  # In any other case, all the parameters are take from `project_constants.jl`
  # file. `input_data` array is initialised to zero and it will not be used.
  else
    include(join([project_name, "_constants.jl"]))
    input_data = zeros(Float64, 1,1)
  end
  # The total simulation time `total_time` is calculated as number of cycles
  # times the cardiac period in seconds.
  const total_time = cycles*cardiac_period

  #Initialise heart data structure
  heart_data =  Heart(inlet_BC_switch,
                            inlet_type,
                            cardiac_period, sys_T,
                            initial_flow, flow_amplitude,
                            input_data)

  inlets = [heart_data]
  if number_of_inlets > 1
      for inlet_number = 2:number_of_inlets
          input_data = loadInputData(join([project_name, "_", inlet_number]))
          push!(inlets, loadGlobalConstants(project_name, inlet_number, inlet_BC_switch :: Int64, inlet_type :: String))
      end
  end

  # Kinematic viscosity `nu` and viscous resitance term `Cf` are calculated as
  # $$
  #   \nu = \frac{\mu}{\rho}, \quad C_f = 8\pi\nu.
  # $$
  # where `Cf` is computed by assuming a parabolic velocity profile
  # along the radial direction.
  const nu = mu/rho
  const Cf = 8*pi*nu
  rho_inv = 1/rho
  viscT = 2*(gamma_profile + 2)*pi*mu

  #Initialise blood data structure
  blood_data =  Blood(mu, rho, rho_inv, Cf, gamma_profile, viscT)
  # --------------------------------------------------------------------------
  # Returns:
  # -------------- -----------------------------------------------------------
  # `heart_data`   `::Heart`
  #
  # `blood_data`   `::Blood`
  #
  # `total_time`   `::Float`
  # --------------------------------------------------------------------------
  return inlets, blood_data, total_time
end

function loadGlobalConstants(project_name, inlet_number,
                             inlet_BC_switch :: Int64,
                             inlet_type :: String)

    input_data = loadInputData(join([project_name, "_", inlet_number]))

    cardiac_period = input_data[end, 1]
    sys_T = 0.0
    initial_flow = 0.0
    flow_amplitude = 0.0

    #Initialise heart data structure
    heart_data = Heart(inlet_BC_switch,
                            inlet_type,
                            cardiac_period, sys_T,
                            initial_flow, flow_amplitude,
                            input_data)

      return heart_data
end
