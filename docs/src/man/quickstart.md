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
        verbose::Bool = false,
        out_files::Bool = false,
        conv_ceil::Bool = false,
        )
```

whose main argument is the name of a `.yaml` file (see [here](https://learnxinyminutes.com/docs/yaml/) for an introduction to yaml) containining the description of the vascular network, blood properties, and numerical solver parameters (see [Configuration](config.md)).

Running a simulation consists in calling

```julia
using openBF
run_simulation("input.yml", verbose=true)
```

This will create a `<project_name>_results` folder containing all the output files from the simulation.
