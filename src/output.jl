flush_to_temp(t::Float64, n::Network, config) = flush_to_temp.(t, values(n.vessels), Ref(config))
function flush_to_temp(t::Float64, v::Vessel, config)
    for l in config["write_results"]
        temp = open(v.label * "_$l.temp", "a")
        if l == "P"
            println(
                temp,
                join((t, pressure(v.A[1], v.A0[1], v.beta[1], v.Pext),
                    pressure(v.A[v.node2], v.A0[v.node2], v.beta[v.node2], v.Pext),
                    pressure(v.A[v.node3], v.A0[v.node3], v.beta[v.node3], v.Pext),
                    pressure(v.A[v.node4], v.A0[v.node4], v.beta[v.node4], v.Pext),
                    pressure(v.A[end], v.A0[end], v.beta[end], v.Pext)), " "),
            )
        elseif l == "Q"
            println(
                temp,
                join((t, v.Q[1], v.Q[v.node2], v.Q[v.node3], v.Q[v.node4], v.Q[end]), " "),
            )
        elseif l == "A"
            println(
                temp,
                join((t, v.A[1], v.A[v.node2], v.A[v.node3], v.A[v.node4], v.A[end]), " "),
            )
        elseif l == "u"
            println(
                temp,
                join((t, v.u[1], v.u[v.node2], v.u[v.node3], v.u[v.node4], v.u[end]), " "),
            )
        end
        close(temp)
    end
end


move_temp_to_last(n::Network, config) = move_temp_to_last.(values(n.vessels), Ref(config))
function move_temp_to_last(v::Vessel, config)
    for l in config["write_results"]
        rm(v.label * "_$l.last", force = true)
        mv(v.label * "_$l.temp", v.label * "_$l.last", force = true)
    end
end

append_last_to_out(n::Network, config) = append_last_to_out.(values(n.vessels), Ref(config))
function append_last_to_out(v::Vessel, config)
    for l in config["write_results"]
        out_file = open(v.label * "_$l.out", "a")
        last_a = readdlm(v.label * "_$l.last")
        writedlm(out_file, last_a, " ")
        close(out_file)
    end
end
