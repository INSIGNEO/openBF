## Installation

The latest Julia binary can be downloaded for all platforms from the official [website](https://julialang.org/downloads/).

openBF can be installed via `Pkg.clone`. Start `julia` in a terminal session and run

```julia
using Pkg
Pkg.clone("https://github.com/INSIGNEO/openBF.git")
```

## Run a simulation

openBF API consists of a single function

```julia
run_simulation(
        yaml_config_path::String;
        verbose::Bool = false,         # show progress bars
        out_files::Bool = false,       # save .out files with the all the cycles
        save_stats::Bool = false,      # save .conv file with simulation stats
        )
```

whose main argument is the name of a `.yaml` file (see [here](https://learnxinyminutes.com/docs/yaml/) for an introduction to yaml) containining the description of the vascular network, blood properties, and numerical solver parameters (see [Configuration](config.md)).

Running a simulation consists in calling

```julia
using openBF
run_simulation("input.yml", verbose=true, save_stats=true)
```

This will create a `<project_name>_results` folder containing all the output files from the simulation.

## Outputs

Simulation results are saved in the `<project_name>_results` directory or in the `output_directory` specified. There, for each artery in your network, you'll find a file per output quantity (these are specified in the `write_results` list).

By default only `.last` files are saved. Thes contains only the _last_ simulated cardiac cycle; at convergence, these are the files you want to look at.

You can also set `out_files = true` when calling `run_simulation`. This results in the writing of `.out` files which contain the _whole_ simulation history, i.e., all the cardiac cycles simulated.

When setting `save_stats=true`, a `<project_name>.conv` file is written in the results dolder. This reads

```
<has the simulation converged> # boolean: true or false
<elapsed cardiac cycles> # integer
<elapsed time> # float, seconds
<memory allocated> # integer, bytes
<garbace collection % time> # float, percentage
```

### `.last`/`.out` Format

These files contain as many rows as defined in the config by the `jump` parameter (default `100`), and 6 columns. The first column contains the simulation time for the current cardiac cycles; column 2-6 report waveforms at five locations along the vessel, namely _inlet_, _1/4th_ of the length, _mid-point_, _3/4th_ of the length and _outlet_.