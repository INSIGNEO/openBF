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
number_of_nodes(config::Vector{Dict{Any,Any}}) = maximum(c["tn"] for c in config)
function Network(
    config::Vector{Dict{Any,Any}},
    blood::Blood,
    heart::Heart,
    Ccfl::Float64;
    verbose = true,
)
    prog = verbose ? Progress(length(config); desc = "Building network:") : nothing

    graph = SimpleDiGraph(number_of_nodes(config))

    vessels = Dict()
    for vessel_config in config
        vessel = Vessel(vessel_config, blood, Ccfl)
        add_edge!(graph, vessel.sn, vessel.tn)
        vessels[(vessel.sn, vessel.tn)] = vessel
        verbose && next!(prog)
    end
    check(graph)
    Network(graph, vessels, blood, heart)
end

function check(g::SimpleDiGraph)
    δin(g) != 0 && error("no input vessel")
    δout(g) != 0 && error("no output vessel(s)")
    has_self_loops(g) && error("self loop detected, i.e. sn == tn")
    Δin(g) > 2 || Δout(g) > 2 && error("vertex belonging to more than 3 vessels")
end
