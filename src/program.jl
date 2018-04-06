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


"""
    runSimulation(input_filename :: String; verbose :: Bool = true, clean :: Bool = false)

Main loop wrapper to start a simulation.
"""
function runSimulation(input_filename :: String;
                       verbose :: Bool = true, out_files :: Bool = true)

    data = loadSimulationFiles(input_filename)
    makeResultsFolder(data)
    blood = buildBlood(data["blood"])

    verbose && println("Build arterial network \n")

    jump = data["solver"]["jump"]
    vessels, edges = buildArterialNetwork(data["network"], blood, jump)
    Ccfl = data["solver"]["Ccfl"]
    heart = vessels[1].heart
    total_time = data["solver"]["cycles"]*heart.cardiac_T
    timepoints = linspace(0, heart.cardiac_T, jump)

    verbose && println("Start simulation \n")

    current_time = 0.0
    passed_cycles = 0

    verbose && (@printf("Solving cardiac cycle no: %d", passed_cycles + 1); tic())

    counter = 1
    while true
        dt = calculateDeltaT(vessels, Ccfl)
        solveModel(vessels, edges, blood, dt, current_time)
        updateGhostCells(vessels)

        if current_time >= timepoints[counter]
            saveTempData(current_time, vessels, counter)
            counter += 1
        end

        if (current_time - heart.cardiac_T*passed_cycles) >= heart.cardiac_T &&
          (current_time - heart.cardiac_T*passed_cycles + dt) > heart.cardiac_T

            # closeTempFiles(vessels)

            if passed_cycles + 1 > 1
                err = checkConvergence(edges, vessels)
                verbose && @printf(" - Error = %4.2f%%\n", err)
            else
                err = 100.0
                verbose && @printf("\n")
            end


            # openCloseLastFiles(vessels)
            transferTempToLast(vessels)

            out_files && transferLastToOut(vessels)

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

    # closeTempFiles(vessels)
    # transferTempToOut(vessels)
    writeResults(vessels)

    # if clean == true
    #     cleanOuts(vessels)
    #     # cleanTemps(vessels)
    # end

    cd("..")
end
