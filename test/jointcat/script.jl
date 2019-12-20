using Revise
using openBF
cd("/home/ivan/repos/openBF/test/single-artery/")
data = openBF.loadSimulationFiles("single-artery.yml")
blood = openBF.buildBlood(data["blood"])
jump = data["solver"]["jump"]
vessels, edges = openBF.buildArterialNetwork(data["network"], blood, jump)
openBF.makeResultsFolder(data)
Ccfl = data["solver"]["Ccfl"]
heart = vessels[1].heart
total_time = data["solver"]["cycles"]*heart.cardiac_T
timepoints = range(0, stop=heart.cardiac_T, length=jump)
current_time = 0.0
passed_cycles = 0
counter = 1

dt = openBF.calculateDeltaT(vessels, Ccfl)
openBF.solveModel(vessels, edges, blood, dt, current_time)
openBF.updateGhostCells(vessels)
current_time += dt
plot(vessels[1].u)
