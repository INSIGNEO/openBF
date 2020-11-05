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

calc_error(w_last, w_tmp) = w_last.-w_tmp
calc_conv(err, w_last) = maximum(abs.(err./w_last))*100.0
calc_norm(err) = sqrt(sum(err.^2))

Base.@kwdef struct ConvError{T}
    max::Float64 = 100.0
    loc::String = ""
    var::String = ""
    prev::T
end

ConvError(::ConvCriteria{:err}) = ConvError(prev=100.0)
ConvError(::ConvCriteria{:norm}, N::Int64) = ConvError(prev=zeros(N,2).+100_000)
ConvError(criteria::String, args...) = ConvError(ConvCriteria(criteria), args...)

"""
    checkConvergence(edge_list, vessels :: Array{Vessel, 1})

Compute the maximum error in the pressure and volumetric flow rate waveforms between two
cardiac cycles at the midpoint of all the vessels in the network. Returns maximum error
and the label of the vessels where this is occurring.
"""
function checkConvergence(::ConvCriteria{:err}, vessels :: Array{Vessel, 1})
    maxerr = -1.0
    maxlbl, maxlblq = "", ""
    @inbounds for v=vessels
        for (w_last, w_temp, lblq) in zip([v.Q_l, v.P_l],[v.Q_t,v.P_t], ["Q","P"])
            err = calc_error(w_last[:,4], w_temp[:,4])
            conv_err = calc_conv(err, w_last[:,4])
            if conv_err > maxerr
                maxerr = conv_err
                maxlbl = v.label
                maxlblq = lblq
            end
        end
    end
    return maxerr, maxlbl, maxlblq, maxerr
end


function calcNorms(vessels :: Array{Vessel, 1})
    norms = zeros(length(vessels),2)
    @inbounds for (i,v) in enumerate(vessels)
        for (w_last, w_temp, col) in zip([v.Q_l, v.P_l], [v.Q_t,v.P_t], [1, 2])
            err = calc_error(w_last[:,4], w_temp[:,4])
            norm = calc_norm(err)
            norms[i,col] = calc_norm(err)
        end
    end
    return norms
end

function checkConvergence(::ConvCriteria{:norm}, vessels :: Array{Vessel, 1},
                              previous_norms::Array{Float64,2})
    current_norms = calcNorms(vessels)
    delta =  current_norms# current_norms.*100 ./previous_norms
    # delta[:,1] ./= maximum(delta[:,1])
    maxnorm, ci = findmax(abs.(delta[:,2]))
    
    maxnormlbl = vessels[ci[1]].label
    maxlblq = ci[1] == 1 ? "Q" : "P"
    return maxnorm, maxnormlbl, maxlblq, current_norms
    # ConvError(max=maxnorm, loc=maxnormlbl, var=maxlblq, prev=current_norms)
end

function checkConvergence(criteria::String, args...)
    # vessels :: Array{Vessel, 1},
    #                       previous_err :: Union{Array{Float64,2},Float64},
    #                       conv_criteria :: String)
    checkConvergence(ConvCriteria(criteria), args...)
end
    # if conv_criteria == "err"
    #     return checkConvergenceError(vessels)
    # elseif conv_criteria == "norm"
    #     return checkConvergenceNorm(vessels, previous_err)
    # end
# end
