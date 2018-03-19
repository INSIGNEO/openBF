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

# The `main.jl` file is where `openBF` is implemented. Only
# [`openBF`](openBF.html) library is needed to be imported.
# Until the official `openBF` Julia `Pkg` is created, the
# library is loaded locally.
# push!(LOAD_PATH, "src/")
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
# project_name = ARGS[1]

parsed_args = openBF.parseCommandline()
project_name = parsed_args["project_name"]
verbose = parsed_args["verbose"]

# In `main.jl` all functions from `openBF` library are called using the dot
# notation: `library.function(parameters)`. Function
# [`projectPreamble`](initialise.html#projectPreamble) checks for all
# the files
# needed by the simulation and warns whether any is missing, it prints
# the initial `openBF` logo, and makes the directory structure to save
# temporary and final results.

# `project_constants.jl` file is user defined and must be in the same folder
# where the simulation is started (see [tutorial](../index.html#tutorial)
# page). Here `project_constants.jl` content is loaded in memory for further
# use.
verbose && println("Loading $project_name files")

# include(join([project_name, "_constants.jl"]))

# openBF.projectPreamble(project_name, no_out, no_inputs, number_of_inlets)

# ### Arterial system
# Arterial model and structure is encoded in a `.csv` file user defined. It
# must be in the same folder in which the simulation is started.
# [`readModelData`](initialise.html#readModelData) reads the
# `.csv` file and fill a 2D matrix with all the informations.
# model = openBF.readModelData(join([project_name, ".csv"]))

# project_inputs = openBF.loadSimulationFiles(project_name)

# Data from `project.csv`, `project_constants.jl`, and `project_inlet.dat` (if
# specified) are used to create instances of [`BTypes`](BTypes.html) data
# structures by [`loadGlobalConstants`](initialise.html#loadGlobalConstants).
# inlets, blood, total_time = openBF.loadGlobalConstants(project_name,
#   inlet_BC_switch, inlet_type, cycles, rho, mu, gamma_profile, number_of_inlets)


inputs = openBF.loadSimulationFiles(project_name)

constants = inputs[1]
model = inputs[2]
inlets = inputs[3]
blood = inputs[4]
total_time = inputs[5]

Ccfl = constants["Ccfl"]
initial_pressure = constants["initial_pressure"]

heart = inlets[1]
# if length(inlets) == 1
#     inlets = inlets[1]
# end

verbose && println("Build arterial network \n")

vessels, edges = openBF.buildArterialNetwork(model, heart, blood)

# Before starting the main loop the counter`current_time` is set to zero. It
# will be updated to keep track of time within the simulation.
verbose && println("Start simulation \n")

current_time = 0.0

# In order to show the progress bar an initial estimate of the total running
# time is needed. Thus, $\Delta t$ is calculated with system initial
# conditions
# and used to compute the number or total iterations before the end of the
# simulation. See [ProgressMeter](https://github.com/timholy/ProgressMeter.jl)
# documentation for `Progress` options.

# The simulation is ran in a `while` loop to be ended whether the simulation
# reached convergence or a user defined finish time.
passed_cycles = 0
verbose && (@printf("Solving cardiac cycle no: %d", passed_cycles + 1); tic())

counter = 1

jump = 100
timepoints = linspace(0, heart.cardiac_T, jump)
while true
  # At the beginning of each time step the $\Delta t$ is computed with
  # [`calculateDeltaT`](godunov.html#calculateDeltaT). This is because
  # during the simulation the local properties of each vessel will change and,
  # as a consequence, also the $\Delta t$ will change. The `current_time` is
  # then calculated by incrementing it with the new `dt`.
  dt = openBF.calculateDeltaT(vessels, Ccfl)

  # [`solveModel`](godunov.html#solveModel) reads the `grafo` object
  # and runs the solver for each part of it. This function can distinguish
  # between inlet, bifurcation, conjunction, anastomosis, and outlet.
  #Solve arteries
  openBF.solveModel(vessels, inlets, edges, blood, dt, current_time)

  # [`updateGhostCells`](boundary_conditions.html#updateGhostCells)
  # updates all vessels ghost cells after the solver ends one iteration.
  openBF.updateGhostCells(vessels)

  # All quantities in each vessel are stored by
  # [`saveTempData`](IOutils.html#saveTempData) in `.temp` files until the
  # end of the cardiac cycle.
  if current_time >= timepoints[counter]
    openBF.saveTempData(current_time, vessels)
    counter += 1
  end

  if (current_time - heart.cardiac_T*passed_cycles) >= heart.cardiac_T &&
      (current_time - heart.cardiac_T*passed_cycles + dt) > heart.cardiac_T

      openBF.closeTempFiles(vessels)

      err = openBF.checkConvergence(edges, vessels, passed_cycles)
      verbose && @printf(" - Error = %4.2f%%\n", err)

      openBF.transferLastToOut(vessels)
      openBF.openCloseLastFiles(vessels)
      openBF.transferTempToLast(vessels)
      openBF.openTempFiles(vessels)

      if err <= 5.
          break
      end

    passed_cycles += 1
    verbose && @printf("Solving cardiac cycle no: %d", passed_cycles + 1)

    timepoints += heart.cardiac_T
    counter = 1
  end

  # When the `current_time` is equal to the `total_time` defined by the user
  # exit from `while` loop. This would also mean that an error smaller than 5%
  # has not been achieved. The main is exited by raising an error containing
  # the error value.
  current_time += dt
  if current_time >= total_time
    verbose && println("Not converged after $passed_cycles cycles, End!")
    break
  end
end
verbose && (@printf "\n"; toc())

# Make sure that data from `.temp` files are transferred.
openBF.closeTempFiles(vessels)
openBF.transferTempToOut(vessels)

clean = parsed_args["clean"]
if clean == true
    cleanOuts(vessels)
    cleanTemps(vessels)
    rm("$project_name.csv")
    rm("$project_name\_inlet.dat")
    rm("$project_name\_constants.yml")
end

cd("..")
# run(`rm main.jl`)
