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

function parse_main_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "input_filename"
            help = ".yml input file name"
            required = true
        "--verbose", "-v"
            help = "Print STDOUT (default false)"
            action = :store_true
        "--out_files", "-f"
            help = "Save complete results story rather than only the last cardiac cycle (default true)"
            action = :store_true
        "--conv_ceil", "-c"
            help = "Ceil convergence value to 100 mmHg (default true)"
            action = :store_true
    end

    return parse_args(s)
end

parsed_args = parse_main_args()

input_filename = parsed_args["input_filename"]
verbose = parsed_args["verbose"]
out_files = parsed_args["out_files"]
conv_ceil = parsed_args["conv_ceil"]

openBF.run_simulation(input_filename, verbose=verbose, out_files=out_files, conv_ceil=conv_ceil)

rm("main.jl")
