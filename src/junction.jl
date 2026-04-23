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

# NOTE: Junction holds preallocated work arrays (x, F, J, dx). It is NOT
# safe to solve the same Junction from multiple threads concurrently.
# If the junction loop is ever parallelised, each thread must own a
# private Junction copy (or these arrays must move to a thread-local
# workspace).

"""
    Junction

A generic n-furcation: k ≥ 2 vessels meeting at a single point with
arbitrary parent (`:outlet`) / daughter (`:inlet`) split.

# Fields
- `id::Int`: unique identifier (the shared node ID) for diagnostics.
- `vessels::Vector{Int}`: indices into the network's vessel list.
- `sides::Vector{Symbol}`: `:outlet` or `:inlet`, per attached vessel.
- `signs::Vector{Int}`: +1 for `:outlet`, −1 for `:inlet`. Precomputed.
- `use_total_pressure::Bool`: pressure-continuity flavour (default `false`
  for static pressure, matching current openBF behaviour).

Work arrays (preallocated to avoid per-step allocation):

- `x::Vector{Float64}`: unknowns, length 2k, ordered [A₁,Q₁,A₂,Q₂,…].
- `F::Vector{Float64}`: residual, length 2k.
- `J::Matrix{Float64}`: Jacobian, 2k × 2k.
- `dx::Vector{Float64}`: Newton update, length 2k.
"""
struct Junction
    id::Int
    vessels::Vector{Int}
    sides::Vector{Symbol}
    signs::Vector{Int}
    use_total_pressure::Bool
    x::Vector{Float64}
    F::Vector{Float64}
    J::Matrix{Float64}
    dx::Vector{Float64}
end

function Junction(id::Int, vessels::Vector{Int}, sides::Vector{Symbol};
                  use_total_pressure::Bool=false)
    length(vessels) == length(sides) || throw(ArgumentError(
        "vessels and sides must have equal length"))
    length(vessels) >= 2 || throw(ArgumentError(
        "a junction needs at least 2 attached vessels"))
    all(s -> s in (:inlet, :outlet), sides) || throw(ArgumentError(
        "sides must be :inlet or :outlet"))
    @assert any(s == :outlet for s in sides) "junction $(id): no outlet vessel"
    @assert any(s == :inlet  for s in sides) "junction $(id): no inlet vessel"
    k = length(vessels)
    signs = [s == :outlet ? +1 : -1 for s in sides]
    Junction(id, vessels, sides, signs, use_total_pressure,
             zeros(2k), zeros(2k), zeros(2k, 2k), zeros(2k))
end
