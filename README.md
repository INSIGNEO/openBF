![openBF](https://alemelis.github.io/openbf.jl/images/openBF.svg)

[![License: LGPL v2.1](https://camo.githubusercontent.com/1ea289570b11d27cbe3dfb3c9fdde5eeabac9751/687474703a2f2f696d672e736869656c64732e696f2f62616467652f6c6963656e73652d4c47504c76322e312d627269676874677265656e2e737667)](http://www.gnu.org/licenses/lgpl-2.1)

**openBF** stable version runs on [Julia v0.3](http://julialang.org/downloads/oldreleases.html) (look also [here](https://caretdashcaret.com/2015/10/13/rolling-back-to-julia-v0-3-11/)) and it requires the following additional packages

- [Graphs.jl](https://github.com/JuliaLang/Graphs.jl)
- [ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl)


To install additional packages, start a Julia interactive session and enter `Pkg.add("package_name")`. For further informations refer to [Julia docs](http://docs.julialang.org/en/release-0.4/manual/packages/).

### Linux installation

* install git
```bash
sudo apt install git
```

* get julia 0.3
```bash
wget "https://julialang.s3.amazonaws.com/bin/linux/x64/0.3/julia-0.3.11-linux-x86_64.tar.gz"
```

* untar and move
```bash
tar xC ~/.local/share/ -f julia-0.3.11-linux-x86_64.tar.gz
mv .local/share/julia-483dbf5279 .local/share/julia
```

* create alias
```bash
echo "alias julia='~/.local/share/julia/bin/julia'" >> ~/.bashrc && source ~/.bashrc
```

* install openBF dependencies
```bash
julia
```
```julia
julia> Pkg.add("Graphs")
julia> Pkg.add("ProgressMeter")
julia> exit()
```

* clone openBF
```bash
git clone https://github.com/INSIGNEO/openBF.git
```

* run test
```bash
cd openbf.jl/test/single-artery
julia main.jl single-artery
```
test output

```julia
Warning: replacing module openBF
Load project single-artery files
Start simulation

Running 100%|██████████████████████████████████████████████████| Time: 0:00:04

elapsed time: 5.016703453 seconds
```

### Tests

To use openBF, clone or download this repository. Go to the chosen test `openBF/test/<testname>` folder and launch the simulation as
```
julia main.jl <testname>
```
This requires the following files to be in the same directory of `main.jl`:

- `<testname>.csv` contains the description of the arterial system;
- `<testname>_veins.csv` contains the description of the venous system;
- `<testname>_constants.jl` contains the simulation global variables;
- `<testname>_inlet.dat` contains the inlet BC.

The file `test/plot.py` is a simple Matplotlib script to print the results.

The `main.jl` file can be copied and used as it is for different simulations. The variables `project_name` must be changed accordingly to the `.csv` files name.

### Docs

For the library documentation check the [Docs](https://INSIGNEO.github.io/openBF/Docs/index.html).
