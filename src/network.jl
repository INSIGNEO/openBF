#=
Copyright 2015-2026 INSIGNEO Institute for in silico Medicine

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
    t::Vector{Float64}
    q::Vector{Float64}
end

load_input_data(project_name::String) = readdlm(project_name * "_inlet.dat")
function Heart(inlet_file::String)
    input_data = readdlm(inlet_file)
    cardiac_period = input_data[end, 1]
    Heart(cardiac_period, input_data, input_data[:,1], input_data[:,2])
end

struct Network
    blood::Blood
    heart::Heart
    Ccfl::Float64
    vessels_vec::Vector{Vessel}
    edge_to_eid::Dict{Tuple{Int,Int}, Int}
    eid_to_edge::Vector{Tuple{Int,Int}}
    parent_eids::Vector{NTuple{2,Int32}}
    child_eids::Vector{NTuple{2,Int32}}
    parent_count::Vector{Int8}
    child_count::Vector{Int8}
    is_inlet::BitVector
    is_outlet::BitVector
    topo_order::Vector{Int32}
    junctions::Vector{Junction}
end
number_of_nodes(config::Vector{Dict{Any,Any}}) = maximum(c["tn"] for c in config)
function build_junctions(vessels::Vector{Vessel})::Vector{Junction}
    by_node = Dict{Int, Vector{Tuple{Int,Symbol}}}()
    for (vid, v) in enumerate(vessels)
        push!(get!(by_node, v.sn, Tuple{Int,Symbol}[]), (vid, :inlet))
        push!(get!(by_node, v.tn, Tuple{Int,Symbol}[]), (vid, :outlet))
    end
    junctions = Junction[]
    for (node, attachments) in by_node
        length(attachments) < 2 && continue
        vids  = [a[1] for a in attachments]
        sides = [a[2] for a in attachments]
        @assert any(s == :outlet for s in sides) "node $node: no outlet vessel"
        @assert any(s == :inlet  for s in sides) "node $node: no inlet vessel"
        k = length(vids)
        push!(junctions, Junction(node, vids, sides; use_total_pressure=(k == 2)))
    end
    sort!(junctions, by=j -> j.id)
    return junctions
end

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
    vessels_vec = Vessel[]
    edge_to_eid = Dict{Tuple{Int,Int}, Int}()
    eid_to_edge = Tuple{Int,Int}[]
    for (eid, vessel_config) in enumerate(config)
        vessel = Vessel(vessel_config, blood, jump, tokeep)
        Graphs.add_edge!(graph, vessel.sn, vessel.tn)
        key = (vessel.sn, vessel.tn)
        vessels[key] = vessel
        push!(vessels_vec, vessel)
        edge_to_eid[key] = eid
        push!(eid_to_edge, key)
        verbose && next!(prog)
    end
    check(graph)
    n_vessels = length(vessels_vec)

    parent_eids  = fill((Int32(0), Int32(0)), n_vessels)
    child_eids   = fill((Int32(0), Int32(0)), n_vessels)
    parent_count = zeros(Int8, n_vessels)
    child_count  = zeros(Int8, n_vessels)

    for eid in 1:n_vessels
        (sn, tn) = eid_to_edge[eid]
        # parents: edges whose dst == sn
        pars = [edge_to_eid[(s, sn)] for s in Graphs.inneighbors(graph, sn)
                if haskey(edge_to_eid, (s, sn))]
        parent_count[eid] = Int8(length(pars))
        if length(pars) >= 1; parent_eids[eid] = (Int32(pars[1]), length(pars) >= 2 ? Int32(pars[2]) : Int32(0)); end

        # children: edges whose src == tn
        chs = [edge_to_eid[(tn, d)] for d in Graphs.outneighbors(graph, tn)
               if haskey(edge_to_eid, (tn, d))]
        child_count[eid] = Int8(length(chs))
        if length(chs) >= 1; child_eids[eid] = (Int32(chs[1]), length(chs) >= 2 ? Int32(chs[2]) : Int32(0)); end
    end

    is_inlet  = BitVector(parent_count .== 0)
    is_outlet = BitVector(child_count  .== 0)

    # Build edge-induced DAG for topological sort:
    # edge A → edge B iff dst(A) == src(B), i.e. child_eids
    edge_dag = SimpleDiGraph(n_vessels)
    for eid in 1:n_vessels
        nc = child_count[eid]
        c1, c2 = child_eids[eid]
        nc >= 1 && Graphs.add_edge!(edge_dag, eid, Int(c1))
        nc >= 2 && Graphs.add_edge!(edge_dag, eid, Int(c2))
    end
    topo_order = Int32.(Graphs.topological_sort_by_dfs(edge_dag))

    junctions = build_junctions(vessels_vec)

    Network(blood, heart, Ccfl,
            vessels_vec, edge_to_eid, eid_to_edge,
            parent_eids, child_eids, parent_count, child_count,
            is_inlet, is_outlet, topo_order, junctions)
end

function check(g::SimpleDiGraph)
    δin(g) != 0 && error("no input vessel")
    δout(g) != 0 && error("no output vessel(s)")
    has_self_loops(g) && error("self loop detected, i.e. sn == tn")
    Δin(g) > 2 || Δout(g) > 2 && error("vertex belonging to more than 3 vessels")
end
