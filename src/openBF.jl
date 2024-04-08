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

__precompile__(true)
module openBF

using Printf
using LinearAlgebra
using DelimitedFiles
using StaticArrays
using YAML
using Glob
using Graphs
using ProgressMeter
using Statistics

export run_simulation

include("vessel.jl")
include("network.jl")
include("solver.jl")
include("simulation.jl")
include("boundary_conditions.jl")
include("conjunctions.jl")
include("bifurcations.jl")
include("output.jl")
include("anastomosis.jl")
end
