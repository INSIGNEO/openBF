
push!(LOAD_PATH, "./src2")
using openBF

parsed_args = openBF.parseCommandline()

input_filename = parsed_args["input_filename"]
verbose = parsed_args["verbose"]
out_files = parsed_args["out_files"]
conv_ceil = parsed_args["conv_ceil"]

openBF.runSimulation(input_filename, verbose=verbose, out_files=out_files, conv_ceil=conv_ceil)

# rm("main.jl")

#-------------------------------------------
openBF.projectPreamble(project_name)

println("Load project $project_name files")
include(join([project_name, "_constants.jl"]))

model = openBF.readModelData(join([project_name, ".csv"]))

heart, blood_prop, total_time = openBF.loadGlobalConstants(project_name,
  inlet_BC_switch, cycles, rho, mu, gamma_profile)

grafo = Graphs.simple_graph(length(model[:,1])+1)

vessels = [openBF.initialiseVessel(model[1,:], 1, heart, blood_prop,Pext,
  initial_pressure, Ccfl)]

Graphs.add_edge!(grafo, vessels[1].sn, vessels[1].tn)

for i in 2:length(model[:,1])


  push!(vessels, openBF.initialiseVessel(model[i,:], i, heart, blood_prop,Pext,
    initial_pressure, Ccfl))

  Graphs.add_edge!(grafo, vessels[end].sn, vessels[end].tn)
end


edgess = Graphs.edges(grafo)


println("Start simulation \n")
current_time = 0


dts  = zeros(Float64, length(edgess))
dt = openBF.calculateDeltaT(vessels, dts)


prog = ProgressMeter.Progress(int(total_time/dt)*2, 1, "Running ", 50)


passed_cycles = 0

tic()
counter = 0
counter_v = 0
while true

  dt = openBF.calculateDeltaT(vessels, dts)

  if venous_model
    dt_v = openBF.calculateDeltaT(vessels_v, dts_v)
    dt = minimum([dt, dt_v])
  end

  current_time += dt


  openBF.solveModel(grafo, edgess, vessels, heart,
                    blood_prop, dt, current_time)

  openBF.updateGhostCells(vessels)

  if counter == 100

    openBF.saveTempData(current_time, vessels)
    counter = 0
  else
    counter += 1
  end

  ProgressMeter.next!(prog)

  if (current_time - heart.cardiac_T*passed_cycles) >= heart.cardiac_T &&
      (current_time - heart.cardiac_T*passed_cycles + dt) > heart.cardiac_T

      openBF.closeTempFiles(vessels)
      openBF.transferLastToOut(vessels)
      openBF.openCloseLastFiles(vessels)
      openBF.transferTempToLast(vessels)
      openBF.openTempFiles(vessels)


    passed_cycles += 1

    openBF.openCloseLastFiles(vessels)
    openBF.transferTempToLast(vessels)
    openBF.openTempFiles(vessels)
  end

  if current_time >= total_time

    break
  end
end
@printf "\n"
toc()


openBF.closeTempFiles(vessels)
openBF.transferTempToOut(vessels)

if venous_model
  openBF.closeTempFiles(vessels_v)
  openBF.transferTempToOut(vessels_v)
end

cd("..")
