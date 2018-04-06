#=
Copyright 2018 INSIGNEO Institute for in silico Medicine

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


function checkConvergence(edge_list, vessels :: Array{Vessel, 1})
    qss = ["_Q", "_P"]

    m_err = []
    for qs in qss
        err = zeros(size(edge_list)[1]).+100

        for i in 1:size(edge_list)[1]
            v = vessels[i]
            lbl = v.label

            filename_last = join([lbl, qs, ".last"])
            filename_temp = join([lbl, qs, ".temp"])
            w_last = readdlm(filename_last)
            w_temp = readdlm(filename_temp)

            # if length(w_last[:,1]) == length(w_temp[:,1])
            #     err[i] = maximum(abs.((w_last[2:end,:].-w_temp[2:end,:])./w_last[2:end,:])*100)
            # end
            err[i] = maximum(abs.((w_last[:,4].-w_temp[:,4])./w_last[:,4])*100)

        end
        push!(m_err, maximum(err))
    end
    return maximum(m_err)
end
