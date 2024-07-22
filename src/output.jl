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

function save_waveforms(idx::Int64, t::Float64, v::Vessel)
    for k in keys(v.waveforms)
        if k == "A"
            v.waveforms[k][idx, :] = [t, v.A[1], v.A[v.node2], v.A[v.node3], v.A[v.node4], v.A[end]]
        elseif k == "Q"
            v.waveforms[k][idx, :] = [t, v.Q[1], v.Q[v.node2], v.Q[v.node3], v.Q[v.node4], v.Q[end]]
        elseif k == "u"
            v.waveforms[k][idx, :] = [t, v.u[1], v.u[v.node2], v.u[v.node3], v.u[v.node4], v.u[end]]
        elseif k == "P"
            # TODO
            # this ---------------------------v should be a call to compute_pressure instead and handle time
            # v.waveforms[k][idx, :] = [t, [pressure(v.A[i], v.A0[i], v.beta[i], v.Pext) - v.Pout for i=(1, v.node2, v.node3, v.node4, v.M)]...]
        end
    end
end


function flush_waveforms(v::Vessel)
    for k in keys(v.waveforms)
        open("$(v.label)_$k.last", "w") do io
            writedlm(io, v.waveforms[k], " ")
        end
    end
end


function append_last_to_out(v::Vessel)
    for k in keys(v.waveforms)
        open("$(v.label)_$k.out", "a") do io
            last = readdlm(v.label * "_$k.last")
            writedlm(io, last, " ")
        end
    end
end
