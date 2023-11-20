struct Heart
    cardiac_period::Float64
    input_data::Array{Float64,2}
end

load_input_data(project_name::String) = readdlm(project_name * "_inlet.dat")

function Heart(project_name::String)
    input_data = load_input_data(project_name)
    cardiac_period = input_data[end, 1]
    Heart(cardiac_period, input_data)
end

struct Network
    graph::SimpleDiGraph
    vessels::Dict{Tuple{Int,Int},Vessel}
    blood::Blood
    heart::Heart
end

function Network(
    config::Vector{Dict{Any,Any}},
    blood::Blood,
    heart::Heart,
    Ccfl::Float64;
    verbose = false,
)
    prog = verbose ? Progress(length(config); desc = "Building network:") : Nothing

    graph = SimpleDiGraph(length(config) + 1)

    vessels = Dict()
    for (i, vessel_config) in enumerate(config)
        vessel = Vessel(vessel_config, i, blood, Ccfl)
        add_edge!(graph, vessel.sn, vessel.tn)
        vessels[(vessel.sn, vessel.tn)] = vessel
        verbose && next!(prog)
    end
    Network(graph, vessels, blood, heart)
end
