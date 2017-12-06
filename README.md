![openBF](https://alemelis.github.io/openbf.jl/images/openBF.svg)

[![INSIGNEO](https://img.shields.io/badge/-INSIGNEO-red.svg?logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAAA4AAAAOCAQAAAC1QeVaAAABGElEQVQY012QvyvEcRzGH12p4%2BRHKW6Rugz%2BAnUWGQyUVYRBFuVMBgmL%2BBdsZzBQymT1YyDJxiaTG3RX576f1%2FubXPE23I8uz7M9z7vn6XlLkiRPhFmOuObK8t7BbpRRA19D3OF1Rt7GG%2B91Oxojwnm1jZAN49F0ZRjHKYQRFVMUcE4%2BOhs5YYUqt2wzqrCJc%2B%2BJZoeiTKkrTrNuiwqPeJiRJG8Pq%2BzQL0k2hduDiPFytyTZMj%2F2FY4lqdKLE4lv3JOSRI5qgAtJKqZwYvGM22RNsAPL1yaELM6LwhbOuf4hXOJ2qDhthrPWapHDKUV9kmyJX5xTJjz52ROynOH8Mle%2FtAVovs9xsPmWoHiQPW4oU%2BCJfRuoqX8d8dI8uuCeiQAAAABJRU5ErkJggg%3D%3D)](https://insigneo.org/)

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

linux/OSX: [![Tests](https://img.shields.io/badge/julia%20v0.3.11-Tests%20pass-brightgreen.svg)](https://julialang.org/downloads/oldreleases.html)
[![Tests](https://img.shields.io/badge/julia%20v0.6.0-Tests%20pass-brightgreen.svg)](https://julialang.org/downloads/)

openBF is an open-source 1D blood flow solver based on MUSCL finite-volume numerical scheme, written in [Julia](https://julialang.org/downloads/) and released under [Apache 2.0](http://www.apache.org/licenses/LICENSE-2.0) free software license.


### Docs

Check the library [website](https://INSIGNEO.github.io/openBF/Docs/index.html) for documentation and tutorials.


### Models

A collection of 1D networks (from the literature) solved by means of openBF can be found in the [openBF-hub](https://github.com/alemelis/openBF-hub) repository


### Julia and openBF installation

- Obtain latest Julia release for your platform [here](https://julialang.org/downloads/); on linux
```bash
$ cd ~/Dowloads
$ wget https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.0-linux-x86_64.tar.gz
$ tar -xzvf julia-0.6.0-linux-x86_64.tar.gz
$ mv julia-0.6.0-linux-x86_64 ~/julia0.6
$ rm julia-0.6.0-linux-x86_64.tar.gz
$ echo "alias julia='~/julia0.6/bin/julia'" > ~/.bashrc
$ source ~/.bashrc
$ julia
```

- clone this repository and copy openBF source files in `.julia/v0.6/` folder; on linux
```bash
$ git clone https://github.com/INSIGNEO/openBF
$ cd openBF
$ mkdir -p ~/.julia/v0.6/openBF
$ cp -r src ~/.julia/v0.6/openBF/
$ cp main.jl ~/.julia/v0.6/openBF/
$ echo 'run(`rm main.jl`)' >> ~/.julia/v0.6/openBF/main.jl
$ mkdir -p ~/.julia/v0.6/BTypes/src
$ cp src/BTypes.jl ~/.julia/v0.6/BTypes/src/
$ echo "alias openBF='cp ~/.julia/v0.6/openBF/main.jl ./main.jl && julia main.jl $1'" >> ~/.bashrc
$ echo "alias openBF-fast='cp ~/.julia/v0.6/openBF/main.jl ./main.jl && julia --check-bounds=no --math-mode=fast main.jl $1'" >> ~/.bashrc
$ source ~/.bashrc
```

### Tests

To use openBF, clone or download this repository. Go to the chosen test `openBF/test/<testname>` folder and launch the simulation as

```
$ julia main.jl
```

This requires the following files to be in the same directory of `main.jl`:

- `<testname>.csv` contains the description of the arterial system;
- `<testname>_veins.csv` contains the description of the venous system;
- `<testname>_constants.jl` contains the simulation global variables;
- `<testname>_inlet.dat` contains the inlet BC.

The file `test/plot.py` (where available) is a simple matplotlib script to print the results.

The `main.jl` file can be copied and used as it is for different simulations. The variable `project_name` must be changed accordingly to the `.csv` files name.
