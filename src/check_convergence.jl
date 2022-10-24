#=
Copyright 2022 INSIGNEO Institute for in silico Medicine

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

function calcNorms(vessels :: Array{Vessel, 1})
    norms = zeros(length(vessels),2)
    for (i,v) in enumerate(vessels)
        err = v.P_l[:,4] .- v.P_t[:,4]
        norms[i] = sqrt(sum(err.^2))
    end
    return norms
end

"""
    computeConvError(edge_list, vessels :: Array{Vessel, 1})
"""
function computeConvError(vessels :: Array{Vessel, 1})
    current_norms = calcNorms(vessels)
    maxnorm, ci = findmax(current_norms)
    maxnormloc = vessels[ci[1]].label
    return maxnorm, maxnormloc
end

function printConvError(err::Float64, loc::String, conv_ceil::Bool)
    err /= 133.332
    if err > 100.0 && conv_ceil
        @printf(" - Error norm > 100.00 mmHg\n")
    else
        @printf(" - Error norm = %6.2f mmHg @ %s\n", err, loc)
    end

end

checkConvergence(err::Float64, toll::Float64) = err/133.332 <= toll
