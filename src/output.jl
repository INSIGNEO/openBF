flush_to_temp(t::Float64, n::Network) = flush_to_temp.(t, values(n.vessels))
function flush_to_temp(t::Float64, v::Vessel)
    # TODO: config to pick arrays
    for (l, a) in zip(("P", "Q"), (v.P, v.Q))
        temp = open(v.label * "_$l.temp", "a")
        line = join((t, a[1], a[v.node2], a[v.node3], a[v.node4], a[end]), " ")
        println(temp, line)
        close(temp)
    end
end


move_temp_to_last(n::Network) = move_temp_to_last.(values(n.vessels))
function move_temp_to_last(v::Vessel)
    for l in ("P", "Q")
        mv(v.label * "_$l.temp", v.label * "_$l.last", force = true)
    end
end

append_last_to_out(n::Network) = append_last_to_out.(values(n.vessels))
function append_last_to_out(v::Vessel)
    for l in ("P", "Q")
        out_file = open(v.label * "_$l.out", "a")
        last_a = readdlm(v.label * "_$l.last")
        writedlm(out_file, last_a, " ")
        close(out_file)
    end
end
