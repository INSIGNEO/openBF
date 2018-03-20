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

# __openBF__ is a computational library written
# in [Julia](julialang.org) for blood flow simulations
# in straight elastic arteries. Navier-Stokes equations are reduced
# to 1D form and solved along the longitudinal direction. Numerical
# solution is achieved via a first-order finite volume solver
# based on the Godunov' method. For an usage example refer to
# [main.jl](main.html) file and its documentation.
__precompile__()

module openBF

# Instances of `BTypes` data structures can be created as
#
#     BTypes.Vessel(args...)
#     BTypes.Heart(args...)
#     BTypes.Blood(args...)
#
# where the `args...` are specified for each data structure in the
# following. Functions in [initialise.jl](initialise.html) generate
# instances of `BTypes` data structures.
export Vessel, Heart, Blood,
    checkInputFiles, checkInputFiles, copyInputFilesToResultsFolder, loadInletData,
    buildHeart, buildHearts, buildBlood, checkConstants, loadConstants, parseModelRow,
    loadSimulationFiles, meshVessel, detectCapillaries, buildArterialNetwork, buildVessel,
    readModelData, setInletBC, inputFromData, inletCompatibility, setOutletBC,
    outletCompatibility, wk3, newtonSolver, updateGhostCells



# The inlet boundary condition is applied to the first vessel in the
# model. A `Heart` instance contains all inlet information.
#
# --------------------------------------------------------------------
# Arguments:
# ----------------- --------------------------------------------------
# `BC_switch`       Integer flag to select inlet boundary condition.
#                   It takes values from 1 to 4 and it is chosen by
#                   the user in the `project_constants.jl` file.
#
# `cardiac_T`       Cardiac period.
#
# `sys_T`           Length of systole.
#
# `initial_flow`    `Q` (ml/s) initial value to be set within all the vessels.
#
# `flow_amplitude`  `Q` maximum amplitude (ml/s) for the user-defined
#                   inlet boundary condition.
#
# `input_data`      Matrix defined as
#                   `[time::Array, flow::Array]`
#                   to be used as inlet boundary condition.
# --------------------------------------------------------------------
# <a name="Heart"></a>
type Heart
  inlet_type :: String
  cardiac_T :: Float64
  input_data :: Array{Float64,2}
  inlet_number :: Int64
end

# Blood mechanical properties.
#
# --------------------------------------------------------------------
# Arguments:
# -------  -----------------------------------------------------------
# `mu`     Dynamic viscosity (poise)
#
# `rho`    Density (kg/m$^3$)
#
# `Cf`     Viscous loss coefficient (see
#          [`loadGlobalConstants`](initialise.html#loadGlobalConstants)
#          and [`source`](godunov.html#source)).
# --------------------------------------------------------------------
# <a name="Blood"></a>
type Blood
  mu  :: Float64
  rho :: Float64
  rho_inv :: Float64
end

# Each vessel in the arterial system is represented by a single
# instance of type `Vessel`.
#
# --------------------------------------------------------------------
# Arguments:
# ---------------- ---------------------------------------------------
# `label`          Vessel clinical name. It is used in
#                  [`initialiseVessel`](initialise.html#initialiseVessel)
#                  to generate output file names.
#
# `ID`             Vessel ID within the `Graph` structure defined
#                  in [main.jl](main.html#grafo). Unique integer.
#
# `sn`             Source node ID within the `Graph` structure. Unique
#                  integer.
#
# `tn`             Terminal node ID within the `Graph` structure. Unique
#                  integer.
#
# `M`              Number of nodes along the vessel. This is user defined
#                  in the project `project_model.csv` file.
#
# `dx`             $\Delta x$ gives the spatial discretisation and it
#                  is computed once in [initialise.jl](initialise.html#dx).
#
# `Ccfl`           Courant-Friedrichs-Lewy (CFL) number. It is used to
#                  compute the $\Delta t$ by
#                  [`calculateDeltaT`](godunov.html#calculateDeltaT)
#                  at each time step. The CFL is always in the
#                  interval $[0,1]$.
#
# `beta`           Elastic constant used to compute the trans-mural
#                  [`pressure`](converter.html#pressure) `P`.
#
# `gamma`          Elastic constant used by
#                  [`waveSpeed`](converter.html#waveSpeed)
#                  to compute `c`.
#
# `R0`             Unstressed lumen radius.
#
# `A0`             Unstressed cross-sectional area.
#
# `A`              Cross-sectional area at each node along the vessel.
#
# `Q`              Flow-rate at each node along the vessel.
#
# `u`              Axial velocity at each node along the vessel.
#
# `c`              Wave speed at each node along the vessel.
#
# `P`              Trans-mural pressure at each node along the vessel.
#
# `W1M0`           Backward Riemann invariant at $x=L$.
#
# `W2M0`           Forward Riemann invariant.
#
# `U00A`           `A` in the ghost cell between 0 and -1. Ghost cells are
#                  used in [godunov.jl](godunov.html) to compute correctly
#                  the solution at boundaries cells.
#
# `U01A`           `A` in the ghost cell between -1 and -2.
#
# `UM1A`           `A` in the ghost cell between `M` and `M+1`.
#
# `UM2A`           `A` in the ghost cell between `M+1` and `M+2`. Ghost
#                  cells are defined for `Q` in the same manner.
#
# `temp_P_name`    Name of temporary file where to store `P` values
#                  while running the solution. Temporary files (`.temp`) are
#                  created for all the variables (`Q`, `A`, `u`, `c`).
#                  Temporary files are handled by
#                  functions defined in [IOutilities.jl](IOutilities.html).
#
# `out_P_name`     Name of final output file where to store `P` values.
#                  Every variable has its own `.out` file.
#
# `temp_P`         Julia handle to direct the output writing operation.
#                  Every variable has its own `IOStream`.
#
# `Rt`             Outlet reflection coefficient. See
#                  [boundary_conditions.jl](boundary_conditions.html).
#
# `R1`             Three element windkessel proximal resistance.
#
# `R2`             Peripheral resistance.
#
# `Cc`             Peripheral compliance.
#
# `Pc`             Pressure through peripheral compliance.
# --------------------------------------------------------------------
# <a name="Vessel"></a>
"""
Vessel type
"""
type Vessel
  label :: String

  #Topological notation
  ID :: Int64
  sn :: Int64
  tn :: Int64
  inlet :: Bool
  heart :: Heart

  #Numerical constants
  M       :: Int64
  dx      :: Float64
  invDx   :: Float64
  halfDx  :: Float64

  #Physical constants
  beta  :: Array{Float64,1}
  gamma :: Array{Float64,1}
  gamma_ghost :: Array{Float64,1}
  half_beta_dA0dx :: Array{Float64,1}
  A0    :: Array{Float64,1}
  inv_A0    :: Array{Float64,1}
  s_inv_A0    :: Array{Float64,1}
  dA0dx :: Array{Float64,1}
  dTaudx :: Array{Float64,1}
  Pext  :: Float64
  viscT :: Float64

  #Iterative solution
  A :: Array{Float64,1}
  Q :: Array{Float64,1}
  u :: Array{Float64,1}
  c :: Array{Float64,1}
  P :: Array{Float64,1}

  #Riemann invariants
  W1M0 :: Float64
  W2M0 :: Float64

  #Ghost cells
  U00A :: Float64
  U00Q :: Float64
  U01A :: Float64
  U01Q :: Float64

  UM1A :: Float64
  UM1Q :: Float64
  UM2A :: Float64
  UM2Q :: Float64

  #Temporary files name
  temp_P_name :: String
  temp_Q_name :: String
  temp_A_name :: String
  temp_c_name :: String
  temp_u_name :: String

  last_P_name :: String
  last_Q_name :: String
  last_A_name :: String
  last_c_name :: String
  last_u_name :: String

  #Output files name
  out_P_name :: String
  out_Q_name :: String
  out_A_name :: String
  out_c_name :: String
  out_u_name :: String

  #Temporary files IOstreams
  temp_P :: IOStream
  temp_Q :: IOStream
  temp_A :: IOStream
  temp_c :: IOStream
  temp_u :: IOStream

  last_P :: IOStream
  last_Q :: IOStream
  last_A :: IOStream
  last_c :: IOStream
  last_u :: IOStream

  #Saving locations
  node2 :: Int64
  node3 :: Int64
  node4 :: Int64

  #Peripheral boundary condition
  Rt :: Float64

  R1 :: Float64
  R2 :: Float64
  Cc :: Float64
  Pc :: Float64

  #Slope
  slope :: Array{Float64,1}

  #MUSCLArrays
  flux :: Array{Float64,2}
  uStar :: Array{Float64,2}

  vA :: Array{Float64,1}
  vQ :: Array{Float64,1}

  dU :: Array{Float64,2}

  slopesA :: Array{Float64,1}
  slopesQ :: Array{Float64,1}

  Al :: Array{Float64,1}
  Ar :: Array{Float64,1}

  Ql :: Array{Float64,1}
  Qr :: Array{Float64,1}

  Fl :: Array{Float64,2}
  Fr :: Array{Float64,2}

  outlet :: String
end

# ### Import openBF' files

# OpenBF' own types are contained in [BTypes.jl](BTypes.html) where
# data structures for vessels, blood, and numerical scheme
# are specified.
# using BTypes
using YAML
using ArgParse

# Data structures and output files are initialised at the beginning
# of each simulation. [initialise.jl](initialise.html) contains all
# the functions needed to read model description and the
# global settings.
include("initialise.jl")

# Inlet and outlets boundary conditions are applied at each time step,
# and vessels extremity nodes are updated with functions defined in
# [boundary_conditions.jl](boundary_conditions.html).
include("boundary_conditions.jl")

# Interface problems (junctions and bifurcations) are solved via the
# method of characteristics. The function in
# [junctions.jl](junctions.html) discriminates between
# [conjunctions](conjunctions.html) and
# [bifurcations](bifurcations.jl), and calls the appropriate solver.
include("junctions.jl")
include("conjunctions.jl")
include("bifurcations.jl")

# [godunov.jl](godunov.html) contains Godunov' method functions.
include("godunov.jl")
include("solver.jl")

# [MUSCL.jl](MUSCL.html) is where the MUSCL method is implemented.
include("MUSCL.jl")

# Conversions from pressure to cross-sectional area (and vice versa),
# or from cross-sectional area to wave speed (and vice versa) are
# handled by functions in [converter.jl](converter.html); Riemann
# invariants are computed in the same location as well.
include("converter.jl")

# Input and output writing functions are stored defined in
# [IOutils.jl](IOutils.html).
include("IOutils.jl")

# Convergence check functions are all defined in [check_convergence.jl](check_convergence.html).
include("check_convergence.jl")

include("anastomosis.jl")



end
