#=
Copyright 2018 INSIGNEO Institute for in silico Medicine

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

using openBF

parsed_args = openBF.parseCommandline()
input_filename = parsed_args["input_filename"]
verbose = parsed_args["verbose"]

verbose && println("Loading $input_filename")

data = openBF.loadSimulationFiles(input_filename)
openBF.makeResultsFolder(data)

blood = openBF.buildBlood(data["blood"])

verbose && println("Build arterial network \n")

vessels, edges = openBF.buildArterialNetwork(data["network"], blood)

Ccfl = data["solver"]["Ccfl"]
heart = vessels[1].heart
total_time = data["solver"]["cycles"]*heart.cardiac_T
jump = data["solver"]["jump"]
timepoints = linspace(0, heart.cardiac_T, jump)

verbose && println("Start simulation \n")

current_time = 0.0
passed_cycles = 0

verbose && (@printf("Solving cardiac cycle no: %d", passed_cycles + 1); tic())

counter = 1
while true
  dt = openBF.calculateDeltaT(vessels, Ccfl)
  openBF.solveModel(vessels, edges, blood, dt, current_time)
  openBF.updateGhostCells(vessels)

  if current_time >= timepoints[counter]
    openBF.saveTempData(current_time, vessels)
    counter += 1
  end

  if (current_time - heart.cardiac_T*passed_cycles) >= heart.cardiac_T &&
      (current_time - heart.cardiac_T*passed_cycles + dt) > heart.cardiac_T

      openBF.closeTempFiles(vessels)

      if passed_cycles + 1 > 1
          err = openBF.checkConvergence(edges, vessels, passed_cycles)
          verbose && @printf(" - Error = %4.2f%%\n", err)
    else
        err = 100.0
        verbose && @printf("\n")
    end

      openBF.transferLastToOut(vessels)
      openBF.openCloseLastFiles(vessels)
      openBF.transferTempToLast(vessels)
      openBF.openTempFiles(vessels)

      if err <= data["solver"]["convergence tollerance"]
          break
      end

    passed_cycles += 1
    verbose && @printf("Solving cardiac cycle no: %d", passed_cycles + 1)

    timepoints += heart.cardiac_T
    counter = 1
  end

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
end

cd("..")
# run(`rm main.jl`)
