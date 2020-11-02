#=
Copyright 2020 INSIGNEO Institute for in silico Medicine

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
    checkConvergence(edge_list, vessels :: Array{Vessel, 1})

Compute the maximum error in the pressure and volumetric flow rate waveforms between two
cardiac cycles at the midpoint of all the vessels in the network. Returns maximum error
and the label of the vessels where this is occurring.
"""
function checkConvergence(edge_list, vessels :: Array{Vessel, 1})
    maxerr, maxnorm = -1.0, -1.0
    maxlbl, maxnormlbl = "", ""
    @inbounds for i in 1:size(edge_list)[1]
        v = vessels[i]
        lbl = v.label
        
        for (w_last, w_temp) in zip([v.Q_l, v.P_l],[v.Q_t,v.P_t])
            err = w_last[:,4].-w_temp[:,4]
            conv_err = maximum(abs.(err./w_last[:,4])*100)
            if conv_err > maxerr
                maxerr = conv_err
                maxlbl = lbl
            end

            norm = sqrt(sum(err.^2))/maximum(abs.(err))
            if norm > maxnorm
                maxnorm = norm
                maxnormlbl = lbl
            end
        end
    end

    return maxerr, maxlbl, maxnorm, maxnormlbl
end

