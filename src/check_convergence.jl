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

struct ConvCriteria{T} end
ConvCriteria(s::AbstractString) = ConvCriteria{Symbol(s)}()

"""
    computeConvError(edge_list, vessels :: Array{Vessel, 1})

Compute the maximum error in the pressure and volumetric flow rate waveforms between two
cardiac cycles at the midpoint of all the vessels in the network. Returns maximum error
and the label of the vessels where this is occurring.
"""
function computeConvError(::ConvCriteria{:err}, vessels :: Array{Vessel, 1})
    maxerr = -1.0
    maxloc = ""
    @inbounds for v=vessels
        for (w_last, w_temp, lblq) in zip([v.Q_l, v.P_l],[v.Q_t,v.P_t], ["Q","P"])
            err = w_last[:,4] .- w_temp[:,4]
            conv_err = maximum(abs.(err./w_last[:,4]))*100.0
            if conv_err > maxerr
                maxerr = conv_err
                maxloc = v.label
            end
        end
    end
    return maxerr, maxloc
end

function calcNorms(vessels :: Array{Vessel, 1})
    norms = zeros(length(vessels),2)
    @inbounds for (i,v) in enumerate(vessels)
        err = v.P_l[:,4] .- v.P_t[:,4]
        norms[i] = sqrt(sum(err.^2))
    end
    return norms
end

function computeConvError(::ConvCriteria{:norm}, vessels :: Array{Vessel, 1})
    current_norms = calcNorms(vessels)
    maxnorm, ci = findmax(current_norms)
    maxnormloc = vessels[ci[1]].label
    return maxnorm, maxnormloc
end

function computeConvError(criteria::String, args...)
    computeConvError(ConvCriteria(criteria), args...)
end

function printConvError(criteria::String, args...)
    printConvError(ConvCriteria(criteria), args...)
end

function printConvError(::ConvCriteria{:norm}, err::Float64, loc::String)
    err /= 133.332
    if err > 100.0
        @printf(" - Error norm > 100.00 mmHg\n")
    else
        @printf(" - Error norm = %6.2f mmHg @ %s\n", err, loc)
    end

end

function printConvError(::ConvCriteria{:err}, err::Float64, loc::String)
    if err > 100.0
        @printf(" - Conv. error > 100.00%%\n")
    else
        @printf(" - Conv. error = %6.2f%% @ %s\n", err, loc)
    end
end

function checkConvergence(criteria::String, args...)
    checkConvergence(ConvCriteria(criteria), args...)
end

checkConvergence(::ConvCriteria{:err}, err::Float64, toll::Float64)=err <= toll
checkConvergence(::ConvCriteria{:norm}, err::Float64, toll::Float64)=err/133.332 <= toll
