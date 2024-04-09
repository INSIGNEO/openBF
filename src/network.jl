#=
Copyright 2015-2024 INSIGNEO Institute for in silico Medicine

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

struct Heart
    cardiac_period::Float64
    input_data::Array{Float64,2}
end

load_input_data(project_name::String) = readdlm(project_name * "_inlet.dat")
function Heart(inlet_file::String)
    input_data = readdlm(inlet_file)        
    cardiac_period = input_data[end, 1]
    Heart(cardiac_period, input_data)
end

struct Network
    graph::SimpleDiGraph
    vessels::Dict{Tuple{Int,Int},Vessel}
    blood::Blood
    heart::Heart
    Ccfl::Float64
    bifU::MVector{6, Float64}
    bifW::MVector{3, Float64}
    bifF::MVector{6, Float64}
    bifJ::MArray{Tuple{6, 6}, Float64, 2, 36}
    # TODO: add vectors for conj and anastomosis
end
number_of_nodes(config::Vector{Dict{Any,Any}}) = maximum(c["tn"] for c in config)
function Network(
    config::Vector{Dict{Any,Any}},
    blood::Blood,
    heart::Heart,
    Ccfl::Float64,
    jump::Int64,
    tokeep::Vector{String};
    verbose = true,
)
    prog = verbose ? Progress(length(config); desc = "Building network:") : nothing

    graph = SimpleDiGraph(number_of_nodes(config))

    vessels = Dict()
    for vessel_config in config
        vessel = Vessel(vessel_config, blood, jump, tokeep)
        add_edge!(graph, vessel.sn, vessel.tn)
        vessels[(vessel.sn, vessel.tn)] = vessel
        verbose && next!(prog)
    end
    check(graph)

    bifU = @MArray zeros(6)
    bifW = @MArray zeros(3)
    bifF = @MArray zeros(6)
    bifJ = @MArray zeros(6,6) 
    bifJ .+= I(6)
    Network(graph, vessels, blood, heart, Ccfl,
        bifU, bifW, bifF, bifJ)
end

function check(g::SimpleDiGraph)
    δin(g) != 0 && error("no input vessel")
    δout(g) != 0 && error("no output vessel(s)")
    has_self_loops(g) && error("self loop detected, i.e. sn == tn")
    Δin(g) > 2 || Δout(g) > 2 && error("vertex belonging to more than 3 vessels")
end
