using Base.Test
using openBF

# projectPreamble
project_name = "test"
number_of_inlets = 2
first_inlet = "$project_name\_inlet.dat"
inlets = [first_inlet]

@test openBF.checkInputFiles(project_name) == nothing

inlets = openBF.checkInletFiles(project_name, number_of_inlets, inlets)
@test length(inlets) == number_of_inlets

@test openBF.copyInputFilesToResultsFolder(project_name, inlets) == nothing
@test isfile("$project_name.csv")
@test isfile("$project_name\_constants.yml")
@test isfile("$project_name\_inlet.dat")
@test isfile("$project_name\_2_inlet.dat")

model = openBF.readModelData("$project_name.csv")
@test typeof(model) == Array{Any,2}

constants = openBF.loadConstants("$project_name\_constants.yml")
@test typeof(constants) == Dict{Any,Any} #constants

@test openBF.checkConstants(constants) == constants

blood = openBF.buildBlood(constants)
@test typeof(blood) == Blood

inlet_data = openBF.loadInletData("$project_name\_inlet.dat")
@test typeof(inlet_data) == Array{Float64, 2}

heart = openBF.buildHeart(constants, inlet_data, 1)
@test typeof(heart) == Heart

hearts = openBF.buildHearts(constants, inlets)
@test typeof(hearts) == Array{Any, 1}

cd("..")

inputs = openBF.loadSimulationFiles(project_name)
@test typeof(inputs[1]) == Dict{Any,Any} #constants
@test typeof(inputs[2]) == Array{Any,2} #model
@test typeof(inputs[3]) == Array{Any,1} #hearts
@test typeof(inputs[4]) == Blood
@test typeof(inputs[5]) == Float64 #total_time

rm("$project_name\_results", recursive=true)
