# Guide

## Installation

The latest Julia binary can be downloaded for all platforms from the official [website](https://julialang.org/downloads/). Additional installation instructions for Windows can be found [here](http://wallyxie.com/weblog/adding-julia-windows-path-command-prompt/).

openBF is not yet a registereed package but,provided you already have a Julia installation, it can be installed via `Pkg.clone`

```julia
Pkg.clone("https://github.com/INSIGNEO/openBF.git")
```

You can also create (MacOSX/Linux only) an openBF alias as

```bash
$ echo "alias openBF='cp ~/.julia/v0.6/openBF/main.jl ./main.jl && julia main.jl $1'" >> ~/.bashrc
$ source ~/.bashrc
```

## Usage

openBF API consists of a single function

```@docs
openBF.runSimulation
```

whose main argument is the name of a `.yml` file (see [here](https://learnxinyminutes.com/docs/yaml/) for an introduction to yaml) containining the description of the vascular network, blood properties, and numerical solver parameters (see [Configuration file](config.md) section for more information on how to write a configuration file).

Running a simulation is as easy as writing

```julia
using openBF
openBF.runSimulation("input.yml", verbose=true)
```

See [Examples](examples.md) section for detailed usage examples.
