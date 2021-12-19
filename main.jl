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


using openBF

parsed_args = openBF.parseCommandline()

input_filename = parsed_args["input_filename"]
verbose = parsed_args["verbose"]
out_files = parsed_args["out_files"]
conv_ceil = parsed_args["conv_ceil"]

openBF.runSimulation(input_filename, verbose=verbose, out_files=out_files, conv_ceil=conv_ceil)

rm("main.jl")
