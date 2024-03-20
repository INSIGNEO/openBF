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

flush_to_temp(t::Float64, n::Network, tosave::Vector{String}) =
    flush_to_temp.(t, values(n.vessels), Ref(tosave))
function flush_to_temp(t::Float64, v::Vessel, tosave::Vector{String})
    for (l, a) in zip(("P", "Q", "A", "u"), (v.P, v.Q, v.A, v.u))
        l ∉ tosave && continue
        temp = open(v.label * "_$l.temp", "a")
        line = join((t, a[1], a[v.node2], a[v.node3], a[v.node4], a[end]), " ")
        println(temp, line)
        close(temp)
    end
end


move_temp_to_last(n::Network, tosave::Vector{String}) =
    move_temp_to_last.(values(n.vessels), Ref(tosave))
function move_temp_to_last(v::Vessel, tosave::Vector{String})
    for l in ("P", "Q", "A", "u")
        l ∉ tosave && continue

        mv(v.label * "_$l.temp", v.label * "_$l.last", force = true)
    end
end

append_last_to_out(n::Network, tosave::Vector{String}) =
    append_last_to_out.(values(n.vessels), Ref(tosave))
function append_last_to_out(v::Vessel, tosave::Vector{String})
    for l in ("P", "Q", "A", "u")
        l ∉ tosave && continue
        out_file = open(v.label * "_$l.out", "a")
        last_a = readdlm(v.label * "_$l.last")
        writedlm(out_file, last_a, " ")
        close(out_file)
    end
end
