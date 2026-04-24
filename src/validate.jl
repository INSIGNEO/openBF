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

"""
    validate_network(yaml_path) -> Bool

Run topology checks on the network described by `yaml_path`.  Prints a
summary of every problem found and returns `true` if the network is clean,
`false` otherwise.

Checks performed:
- Self-loops (sn == tn for any vessel).
- No inlet vessel (every vessel has at least one parent).
- No outlet vessel (every vessel has at least one child).
- Junctions with all-outlet or all-inlet attachments (no flow through the node).
- Suspiciously high-k junctions (k > 3; warn only, not an error).
- Disconnected subgraphs (network is not a single weakly-connected component).
"""
function validate_network(yaml_path::String)::Bool
    config = YAML.load_file(yaml_path)
    network_cfg = config["network"]

    errors  = String[]
    warnings = String[]

    # ── 1. build a simple adjacency picture ──────────────────────────────────
    by_node = Dict{Int, Vector{Tuple{String,Symbol}}}()   # node → [(label, side)]
    for v in network_cfg
        label = get(v, "label", "?")
        sn    = v["sn"]
        tn    = v["tn"]

        if sn == tn
            push!(errors, "vessel '$label': self-loop (sn == tn == $sn)")
            continue
        end

        push!(get!(by_node, sn, Tuple{String,Symbol}[]), (label, :inlet))
        push!(get!(by_node, tn, Tuple{String,Symbol}[]), (label, :outlet))
    end

    # ── 2. inlet / outlet existence ──────────────────────────────────────────
    has_network_inlet  = false
    has_network_outlet = false
    for (_, atts) in by_node
        sides = [a[2] for a in atts]
        all(s == :outlet for s in sides) && (has_network_inlet = true)   # node with only outflows = source
        all(s == :inlet  for s in sides) && (has_network_outlet = true)  # node with only inflows  = sink
    end
    has_network_inlet  || push!(errors,  "no inlet vessel found (network has no source node)")
    has_network_outlet || push!(errors,  "no outlet vessel found (network has no sink node)")

    # ── 3. per-junction checks ───────────────────────────────────────────────
    for (node, atts) in by_node
        length(atts) < 2 && continue          # leaf nodes are not junctions
        sides  = [a[2] for a in atts]
        labels = [a[1] for a in atts]
        k      = length(atts)

        if all(s == :outlet for s in sides)
            push!(errors, "node $node: all-outlet junction (vessels $(join(labels, ", ")))")
        elseif all(s == :inlet for s in sides)
            push!(errors, "node $node: all-inlet junction (vessels $(join(labels, ", ")))")
        end

        if k > 3
            push!(warnings, "node $node: $k-vessel junction — verify this is intentional")
        end
    end

    # ── 4. connectivity (weak) ───────────────────────────────────────────────
    n_vessels = length(network_cfg)
    if n_vessels > 0
        adj = [Set{Int}() for _ in 1:n_vessels]
        for (i, v) in enumerate(network_cfg)
            sn, tn = v["sn"], v["tn"]
            for (j, w) in enumerate(network_cfg)
                i == j && continue
                if w["sn"] == tn || w["tn"] == tn || w["sn"] == sn || w["tn"] == sn
                    push!(adj[i], j)
                end
            end
        end
        visited = falses(n_vessels)
        stack   = [1]
        while !isempty(stack)
            node = pop!(stack)
            visited[node] && continue
            visited[node] = true
            for nb in adj[node]
                visited[nb] || push!(stack, nb)
            end
        end
        if !all(visited)
            n_unreachable = count(!, visited)
            push!(errors, "$n_unreachable vessel(s) unreachable from vessel 1 (disconnected network)")
        end
    end

    # ── 5. report ────────────────────────────────────────────────────────────
    for w in warnings
        println("  WARNING  $w")
    end
    for e in errors
        println("  ERROR    $e")
    end

    ok = isempty(errors)
    if ok && isempty(warnings)
        println("  OK  network '$(get(config, "project_name", yaml_path))' passed all checks")
    elseif ok
        println("  OK  $(length(warnings)) warning(s), no errors")
    else
        println("  FAIL  $(length(errors)) error(s) found")
    end
    return ok
end
