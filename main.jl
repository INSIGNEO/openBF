#= Copyright (C) 2017 Alessandro Melis.

  This file is part of openBF.

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with open.  If not, see <http://www.gnu.org/licenses/>.
=#

# The `main.jl` file is where `openBF` is implemented. Only
# [`openBF`](openBF.html) library is needed to be imported.
# Until the official `openBF` Julia `Pkg` is created, the
# library is loaded locally.
push!(LOAD_PATH, "src/")
using openBF
# reload("openBF")

# The project name must be specified by the user when launching
# the simulation.
# `openBF` takes care of all the pre-processing and solution steps. A
# simulation
# is started with the following command
#
#     julia main.jl project_name
#
# `project_name` is a string used to initialise all the output files.
project_name = ARGS[1]

# In `main.jl` all functions from `openBF` library are called using the dot
# notation: `library.function(parameters)`. Function
# [`projectPreamble`](initialise.html#projectPreamble) checks for all
# the files
# needed by the simulation and warns whether any is missing, it prints
# the initial `openBF` logo, and makes the directory structure to save
# temporary and final results.
openBF.projectPreamble(project_name)

# `project_constants.jl` file is user defined and must be in the same folder
# where the simulation is started (see [tutorial](../index.html#tutorial)
# page). Here `project_constants.jl` content is loaded in memory for further
# use.
println("Load project $project_name files")
include(join([project_name, "_constants.jl"]))

# ### Arterial system
# Arterial model and structure is encoded in a `.csv` file user defined. It
# must be in the same folder in which the simulation is started.
# [`readModelData`](initialise.html#readModelData) reads the
# `.csv` file and fill a 2D matrix with all the informations.
model = openBF.readModelData(join([project_name, ".csv"]))

# Data from `project.csv`, `project_constants.jl`, and `project_inlet.dat` (if
# specified) are used to create instances of [`BTypes`](BTypes.html) data
# structures by [`loadGlobalConstants`](initialise.html#loadGlobalConstants).
heart, blood_prop, total_time = openBF.loadGlobalConstants(project_name,
  inlet_BC_switch, inlet_type, cycles, rho, mu, gamma_profile)

# The arterial tree is represented as a graph by means of `Graphs` library (
# see [grafo.jl](grafo.html) for a detailed example).
# `simple_graph` creates an empty `GenericGraph` with a user defined number of
# nodes.
#
# ----------------------------------------------------------------------------
# Parameters
# ----------- ----------------------------------------------------------------
# `nodes`     `::Int64` number of total nodes in the graph. Here an estimate
#             of this number is taken by counting `model` matrix rows (one
#             for each vessel) and adding one extra node.
# ----------------------------------------------------------------------------
# Returns
# ----------- ----------------------------------------------------------------
# `graph`     `::GenericGraph` empty graph data structure.
# ----------------------------------------------------------------------------
#
# <a name="grafo"></a>
# grafo = LightGraphs.DiGraph(length(model[:,1])+1)

# Two data structures are used to describe the arterial system. One is a
# collection of [`BTypes`](BTypes.html)`.Vessel` instances, and the second
# is the `grafo` structure. `vessels` collection is filled by using
# [`initialiseVessel`](initialise.html#initialiseVessel). The first
# vessel is inserted by hand.
vessels = [openBF.initialiseVessel(model[1,:], 1, heart, blood_prop,
  initial_pressure, Ccfl)]
edge_list = zeros(Int8, length(model[:,1]), 3)
edge_list[1,1] = vessels[1].ID
edge_list[1,2] = vessels[1].sn
edge_list[1,3] = vessels[1].tn

# The graph is built with `Graphs.add_edge!` function.
#
# ----------------------------------------------------------------------------
# Parameters
# ----------- ----------------------------------------------------------------
# `grafo`     `::GenericGraph` arterial system graph structure.
#
# `sn`        `::Int64` edge source node stored in `Vessel` structure.
#
# `tn`        `::Int64` edge terminal node.
# ----------------------------------------------------------------------------
# LightGraphs.add_edge!(grafo, vessels[1].sn, vessels[1].tn)

# The model matrix is read iteratively starting from the second row, the first
# row contains column headers.
for i in 2:length(model[:,1])

  # Each new `vessel` is instantiated and `push!`ed at the bottom of
  # `vessels` collection.
  push!(vessels, openBF.initialiseVessel(model[i,:], i, heart, blood_prop,
    initial_pressure, Ccfl))
  edge_list[i,1] = vessels[i].ID
  edge_list[i,2] = vessels[i].sn
  edge_list[i,3] = vessels[i].tn
  # LightGraphs.add_edge!(grafo, vessels[end].sn, vessels[end].tn)
end

# # `edge_list` is a list of all the edges in `grafo`
# edge_list = collect(LightGraphs.edges(grafo))
# edge_map = Dict{LightGraphs.SimpleGraphs.SimpleEdge, Int}(e=>i
#   for (i, e) in enumerate(edge_list))
# node_map = Dict{Tuple{Int, Int}, Int}()
# n_i = 1
# for e in edge_list
#   s = LightGraphs.src(e)
#   t = LightGraphs.dst(e)
#   node_map[(s,t)] = n_i
#   n_i += 1
# end

# Before starting the main loop the counter`current_time` is set to zero. It
# will be updated to keep track of time within the simulation.
println("Start simulation \n")
current_time = 0

# In order to show the progress bar an initial estimate of the total running
# time is needed. Thus, $\Delta t$ is calculated with system initial
# conditions
# and used to compute the number or total iterations before the end of the
# simulation. See [ProgressMeter](https://github.com/timholy/ProgressMeter.jl)
# documentation for `Progress` options.
dts  = zeros(Float64, length(edge_list[:,1]))
dt = openBF.calculateDeltaT(vessels, dts)

# # ### Venous system
# # Same as the arterial system
# v_model = join([project_name, "_veins.csv"])
# if (isfile(v_model)) == true
#   model_v = openBF.readModelData(join([project_name, "_veins.csv"]))
#
#   grafo_v = Graphs.simple_graph(length(model_v[:,1])+1)
#
#   vessels_v = [openBF.initialiseVessel(model_v[1,:], 1, heart, blood_prop,
#     initial_pressure, Ccfl)]
#
#   Graphs.add_edge!(grafo_v, vessels_v[1].sn, vessels_v[1].tn)
#
#   for i in 2:length(model_v[:,1])
#
#     push!(vessels_v, openBF.initialiseVessel(model_v[i,:], i, heart, blood_prop,
#       initial_pressure, Ccfl))
#
#     Graphs.add_edge!(grafo_v, vessels_v[end].sn, vessels_v[end].tn)
#   end
#
#   edge_list_v = Graphs.edges(grafo_v)
#
#   dts_v  = zeros(Float64, length(edge_list_v))
#   dt_v = openBF.calculateDeltaT(vessels_v, dts_v)
#
#   dt = minimum([dt, dt_v])
#
#   venous_model = true
# else
#   venous_model = false
# end

prog = ProgressMeter.Progress(Int(ceil(total_time/dt)), 1, "Running ", 50)

# The simulation is ran in a `while` loop to be ended whether the simulation
# reached convergence or a user defined finish time.
passed_cycles = 0

tic()
counter = 0
counter_v = 0
while true
  # At the beginning of each time step the $\Delta t$ is computed with
  # [`calculateDeltaT`](godunov.html#calculateDeltaT). This is because
  # during the simulation the local properties of each vessel will change and,
  # as a consequence, also the $\Delta t$ will change. The `current_time` is
  # then calculated by incrementing it with the new `dt`.
  dt = openBF.calculateDeltaT(vessels, dts)

  # if venous_model
  #   dt_v = openBF.calculateDeltaT(vessels_v, dts_v)
  #   dt = minimum([dt, dt_v])
  # end

  current_time += dt

  # [`solveModel`](godunov.html#solveModel) reads the `grafo` object
  # and runs the solver for each part of it. This function can distinguish
  # between inlet, bifurcation, conjunction, anastomosis, and outlet.
  #Solve arteries
  openBF.solveModel(vessels, heart, edge_list, blood_prop, dt, current_time)

  # [`updateGhostCells`](boundary_conditions.html#updateGhostCells)
  # updates all vessels ghost cells after the solver ends one iteration.
  openBF.updateGhostCells(vessels)

  # All quantities in each vessel are stored by
  # [`saveTempData`](IOutils.html#saveTempData) in `.temp` files until the
  # end of the cardiac cycle.
  if counter == 100
    openBF.saveTempData(current_time, vessels)
    counter = 0
  else
    counter += 1
  end

  # #Solve veins
  # if venous_model
  #   openBF.solveModel(grafo_v, edge_list_v, vessels_v,
  #                     grafo, edge_list, vessels,
  #                     blood_prop, dt, current_time)
  #   openBF.updateGhostCells(vessels_v)
  #
  #   if counter_v == 100
  #     openBF.saveTempData(current_time, vessels_v)
  #     counter_v = 0
  #   else
  #     counter_v += 1
  #   end
  #
  # end

  #Progress bar update
  ProgressMeter.next!(prog)

  # Every time a cardiac cycle has been simulated, this condition returns
  # `true` and data from `.temp` files are transferred to `.out` files (
  # see [IOutils.jl](IOutils.html)).

  if (current_time - heart.cardiac_T*passed_cycles) >= heart.cardiac_T &&
      (current_time - heart.cardiac_T*passed_cycles + dt) > heart.cardiac_T

      openBF.closeTempFiles(vessels)
      openBF.transferLastToOut(vessels)
      openBF.openCloseLastFiles(vessels)
      openBF.transferTempToLast(vessels)
      openBF.openTempFiles(vessels)

      # if venous_model
      #   openBF.closeTempFiles(vessels_v)
      #   openBF.transferLastToOut(vessels_v)
      #   openBF.openCloseLastFiles(vessels_v)
      #   openBF.transferTempToLast(vessels_v)
      #   openBF.openTempFiles(vessels_v)
      # end

    passed_cycles += 1
    # # When at least 3 cardiac cycles have been simulated, waveforms are
    # # checked for convergence.
    # if passed_cycles >= 3

    #   # The error is computed for all vessels in the system by
    #   # [`checkAllQuantities`](check_convergence.html#checkAllQuantities)
    #   # function.
    #   # <a name="check_convergence"></a>
    #   # err = openBF.checkAllQuantities(vessels, passed_cycles, 1000)
    #   err = openBF.checkConvergence(vessels, 5.)

    #   # The convergence is reached when the difference (the error)
    #   # between two consecutive waveforms is less than 5%. In this case, the
    #   # main loop is exited.
    #   if err < 5.
    #     println("\nConverged in $passed_cycles cycles, end!")
    #     break
    #   end
    # end

    # openBF.openCloseLastFiles(vessels)
    # openBF.transferTempToLast(vessels)
    # openBF.openTempFiles(vessels)
  end

  # When the `current_time` is equal to the `total_time` defined by the user
  # exit from `while` loop. This would also mean that an error smaller than 5%
  # has not been achieved. The main is exited by raising an error containing
  # the error value.
  if current_time >= total_time
    # erlog = open("error.log", "w")
    # write(erlog, "Not converged after $passed_cycles cycles, End!")
    # close(erlog)
    break
  end
end
@printf "\n"
toc()

# Make sure that data from `.temp` files are transferred.
openBF.closeTempFiles(vessels)
openBF.transferTempToOut(vessels)

# if venous_model
#   openBF.closeTempFiles(vessels_v)
#   openBF.transferTempToOut(vessels_v)
# end

cd("..")
