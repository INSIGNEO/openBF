# Examples

The `models/` directory contains ready-to-run configurations from published benchmarks.

## Boileau 2015 benchmark models

From _Boileau et al., Int. J. Numer. Methods Biomed. Eng., 2015_:

| Model | Path | Description |
|-------|------|-------------|
| `cca` | `models/boileau2015/cca/` | Single common carotid artery |
| `uta` | `models/boileau2015/uta/` | Upper thoracic aorta |
| `ibif` | `models/boileau2015/ibif/` | Iliac bifurcation |
| `adan56` | `models/boileau2015/adan56/` | 56-artery systemic network |

## Circle of Willis

From _Alastruey et al., J. Biomechanics, 2007_:

| Model | Path | Description |
|-------|------|-------------|
| `circle_of_willis` | `models/alastruey2007/` | 33-vessel cerebral circulation with anastomoses |

## Running an example

```julia
using openBF

cd("models/boileau2015/adan56")
openBF.run_simulation("adan56.yaml")
```

Results appear in `adan56_results/`. Each vessel produces one `.last` file per quantity listed in `write_results`.

## Reading results

```julia
using DelimitedFiles

data = readdlm("adan56_results/aorta_P.last")
# data[:, 1]  — time [s]
# data[:, 2]  — pressure at inlet [Pa]
# data[:, 6]  — pressure at outlet [Pa]
```

## Minimal network

A single straight vessel with a WK3 outlet:

```yaml
project_name: single_vessel
write_results: ["P", "Q"]

solver:
  Ccfl: 0.9
  cycles: 20
  convergence_tolerance: 1.0
  jump: 100

blood:
  rho: 1060.0
  mu: 0.004

network:
  - label: aorta
    sn: 1
    tn: 2
    L: 0.4
    E: 400000.0
    R0: 0.015
    M: 40
    R1: 1.0e7
    R2: 9.0e8
    Cc: 1.0e-9
```
