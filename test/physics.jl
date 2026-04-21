using openBF
using YAML
using DelimitedFiles
using Statistics
using Printf

include(joinpath(@__DIR__, "reference.jl"))

"""
    mass_balance(data, config) -> Dict{node_id => relative_imbalance}

For each internal node (has both incoming and outgoing vessels), compute
|mean(ΣQ_in) - mean(ΣQ_out)| / (0.5 * mean(|ΣQ_in| + |ΣQ_out|)).
Requires Q waveforms in `data`. Returns empty dict if Q is missing.
"""
function mass_balance(data, config)
    vessels = config["network"]

    node_in  = Dict{Int, Vector{String}}()
    node_out = Dict{Int, Vector{String}}()
    for v in vessels
        sn, tn, label = v["sn"], v["tn"], v["label"]
        push!(get!(node_out, sn, String[]), label)
        push!(get!(node_in,  tn, String[]), label)
    end

    result = Dict{Int, Float64}()
    for n in sort(collect(union(keys(node_in), keys(node_out))))
        ins  = get(node_in,  n, String[])
        outs = get(node_out, n, String[])
        isempty(ins) || isempty(outs) && continue  # skip source/sink nodes

        q_in_mats  = [data[l]["Q"][:, 6] for l in ins  if haskey(data, l) && haskey(data[l], "Q")]
        q_out_mats = [data[l]["Q"][:, 2] for l in outs if haskey(data, l) && haskey(data[l], "Q")]
        (isempty(q_in_mats) || isempty(q_out_mats)) && continue

        q_in  = sum(q_in_mats)
        q_out = sum(q_out_mats)
        imbalance = mean(abs.(q_in .- q_out))
        scale = 0.5 * mean(abs.(q_in) .+ abs.(q_out))
        result[n] = scale > 0 ? imbalance / scale : 0.0
    end
    result
end

"""
    flow_continuity(data, config) -> Float64

Ratio of mean outlet flow to mean inlet flow over the last cycle.
Should be ~1.0 (exact only if no compliance; WK3 outlets store flow in Pc).
"""
function flow_continuity(data, config)
    vessels = config["network"]
    sn_set = Set(v["sn"] for v in vessels)
    tn_set = Set(v["tn"] for v in vessels)
    source_nodes = setdiff(sn_set, tn_set)
    sink_nodes   = setdiff(tn_set, sn_set)

    inlet_labels  = [v["label"] for v in vessels if v["sn"] ∈ source_nodes]
    outlet_labels = [v["label"] for v in vessels if v["tn"] ∈ sink_nodes]

    q_in = sum(
        mean(data[l]["Q"][:, 2])
        for l in inlet_labels if haskey(data, l) && haskey(data[l], "Q");
        init=0.0
    )
    q_out = sum(
        mean(data[l]["Q"][:, 6])
        for l in outlet_labels if haskey(data, l) && haskey(data[l], "Q");
        init=0.0
    )
    q_in ≈ 0 ? NaN : q_out / q_in
end

"""
    physics_report(model_name) -> (data, config)

Run the current solver on `model_name`, print mass balance per internal node,
flow continuity, and return the captured data for further use (e.g., waveform_diff).
"""
function physics_report(model_name::String)
    println("=== Physics report: $model_name ===")
    data, config = capture_waveforms(model_name)

    mb = mass_balance(data, config)
    if isempty(mb)
        println("  Mass balance: N/A (no Q waveforms for internal nodes)")
    else
        println("  Mass balance (relative Q imbalance per internal node):")
        for n in sort(collect(keys(mb)))
            @printf "    node %d: %6.4f%%\n" n mb[n]*100
        end
        @printf "    MAX: %6.4f%%\n" maximum(values(mb))*100
    end

    fc = flow_continuity(data, config)
    if isnan(fc)
        println("  Flow continuity: N/A")
    else
        @printf "  Flow continuity (Q_out / Q_in): %.6f\n" fc
    end

    println()
    data, config
end

"""
    physics_check(model_name; mb_tol=0.01, fc_tol=0.01) -> NamedTuple

Returns (mass_balance_ok, flow_continuity_ok, mb_max, fc).
"""
function physics_check(model_name::String; mb_tol=0.01, fc_tol=0.01)
    data, config = capture_waveforms(model_name)
    mb = mass_balance(data, config)
    fc = flow_continuity(data, config)
    mb_max = isempty(mb) ? NaN : maximum(values(mb))
    (
        mass_balance_ok    = isnan(mb_max) || mb_max < mb_tol,
        flow_continuity_ok = isnan(fc)     || abs(fc - 1.0) < fc_tol,
        mb_max             = mb_max,
        fc                 = fc,
    )
end
