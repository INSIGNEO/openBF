# Quickstart

## Installation

Install [Julia](https://julialang.org/downloads/), then add openBF from the REPL:

```julia
julia> ]
(@v1.11) pkg> add https://github.com/INSIGNEO/openBF.git
```

## Running a simulation

The entire API is a single function:

```julia
openBF.run_simulation(
    yaml_config::String;
    verbose::Bool   = true,   # progress bars and convergence output
    out_files::Bool = false,  # write full per-cycle history (.out files)
    save_stats::Bool = false, # write convergence stats (.conv file)
    savedir::String = "",     # output directory (overrides config)
)
```

The only required argument is the path to a `.yaml` configuration file (see [Configuration](config.md)). A minimal call:

```julia
using openBF
openBF.run_simulation("my_network.yaml")
```

Results are written to `<project_name>_results/` in the working directory, or to `savedir` if provided.

Alternatively, use the included `run.jl` script from the terminal:

```bash
julia run.jl my_network.yaml
```

## Validating a network

Before running a simulation you can check the network topology for common errors (self-loops, disconnected subgraphs, malformed junctions) without executing the solver:

```bash
julia run.jl --validate my_network.yaml
```

Exits with status 0 if the network is well-formed, status 1 if any errors are found. Warnings (e.g. unusually high-k junctions) are printed but do not affect the exit code.

## Output files

For each vessel listed in `write_results`, two files are produced:

| File | Contents |
|------|----------|
| `<label>_<field>.last` | Last simulated cardiac cycle only |
| `<label>_<field>.out` | All cycles (only with `out_files=true`) |

Each file has `jump` rows (default 100) and 6 columns:

| Column | Description |
|--------|-------------|
| 1 | Time within the cardiac cycle (s) |
| 2 | Inlet value |
| 3 | Value at ¼ vessel length |
| 4 | Value at mid-point |
| 5 | Value at ¾ vessel length |
| 6 | Outlet value |

Fields: `P` pressure (Pa), `Q` flow (m³/s), `u` velocity (m/s), `A` cross-sectional area (m²).

## Convergence stats

With `save_stats=true`, a `<project_name>.conv` file is written:

```
<converged>        # true / false
<cardiac cycles>   # integer
<elapsed time>     # seconds
<memory allocated> # bytes
<GC time %>        # percentage
```
