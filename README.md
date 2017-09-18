![openBF](https://alemelis.github.io/openbf.jl/images/openBF.svg)

[![License: LGPL v2.1](https://img.shields.io/badge/license-LGPL%202.1-yellow.svg)](http://www.gnu.org/licenses/lgpl-2.1)

linux/OS X: [![Tests](https://img.shields.io/badge/julia%20v0.3.11-Tests%20pass-brightgreen.svg)](https://julialang.org/downloads/oldreleases.html)
[![Tests](https://img.shields.io/badge/julia%20v0.6.0-Tests%20pass-brightgreen.svg)](https://julialang.org/downloads/)


**openBF** runs on [Julia](https://julialang.org/downloads/) and it requires the following additional packages

- [ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl)

To install additional packages, start a Julia interactive session and enter `Pkg.add("package_name")`. For further informations refer to [Julia docs](http://docs.julialang.org/en/release-0.4/manual/packages/).

<!-- #### Install julia v0.6.0-1 with ananaconda
```bash
$ conda config --add channels defaults
$ conda config --add channels conda-forge
$ conda config --add channels brown-data-science
$ conda create -n julia
$ source activate julia
$ conda install -c brown-data-science julia
$ source deactivate
``` -->

### Docs

Check the library [website](https://INSIGNEO.github.io/openBF/Docs/index.html) for documentation and tutorials.

### Tests

To use openBF, clone or download this repository. Go to the chosen test `openBF/test/<testname>` folder and launch the simulation as
```
julia main.jl
```
This requires the following files to be in the same directory of `main.jl`:

- `<testname>.csv` contains the description of the arterial system;
- `<testname>_veins.csv` contains the description of the venous system;
- `<testname>_constants.jl` contains the simulation global variables;
- `<testname>_inlet.dat` contains the inlet BC.

The file `test/plot.py` (where available) is a simple matplotlib script to print the results.

The `main.jl` file can be copied and used as it is for different simulations. The variable `project_name` must be changed accordingly to the `.csv` files name.
