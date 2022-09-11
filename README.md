# __openBF.jl__

[![TUoS](https://img.shields.io/badge/-The%20University%20of%20Sheffield-blue.svg?colorA=ffffff&colorB=009fe3&logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAAA4AAAAOCAMAAAAolt3jAAABsFBMVEUAAABmZplVVYBLPHheUYZJPXlSR3pCOnNgWIdoYI9FPnVMQXc0JGM0KmhEOXI0MW8xMW2Ph6yYkrOWj6%2BOiayinLqWkbIcWJkdWJezrsW2sciclrWmoL0cWpkbXJuCfqQVba4Vba2tqsKBgKWXkbGgm7ioo7%2Bln7umobyppcClob0RdrgSd7ikoLuXkbEeeKwieqmYkrKkn7sAneAAn%2BMBnN8Dn%2BIEj9IEkNQEoOMFi88GicwHnuAInuAVaqkVgrEWa6sWiskYj84ek88mk8wmqOInqOIqhrA4qd45qd87kaA8q988q%2BA9qN0%2Fk51Aqd1CqdxCqd1JnqtNnc9OnKVPnItSmsZUmqZVlsFVnYZVns1bst9ds99ftN9isNtitN9joa1msdtnveZqoaRqrNZqvuZsrdZvk2lylWR0nKZ2oJx%2BfKOBf6WEn5eGr3iIr3eduNWfmbifmrifudWim6mjnKmloLymrz2psTuuqb6uqsOwq8C0taG2tqS4tqy4tra6uTC6ui6%2FsCK%2FsCPBvdHDv9LFwtfGwtjOy9rPzdvUlaHmkpb%2B%2FPz%2B%2Fv7%2F%2FPz%2F%2F%2F9OZpcfAAAAM3RSTlMABQYRExUZHyAgJS8xMTFTVFleZmiJkpOTmp2oqba9wsTFxcfO0dPX3OLj5%2Bf19v39%2Fv4kncL4AAAAvElEQVR4AWNgYBIQV9Q20NdRkOBnZGBg4DExTQ9taIzMMDbhA3PdzMqbmyosXaHc4kCvvHy%2FgAIINzElvL%2B%2FvzcoLQ7E5bX1TOhq624P8bDhBnLZTEwye%2Fr7%2BpNNTFgYgMApJreuprY%2BK94RyOGSlEv19XZx9vFPkpViZ1BRZ7WKdq%2Bsto%2B1ZtZQYhAp1FS2iCopDTaX1yoSZmDgkNHrsAuLcOjUleZkAAGxqtbsnJYyUQYoEFIzMlQVBLEA%2FZgsl9iPrB4AAAAASUVORK5CYII%3D)](https://www.sheffield.ac.uk)
[![INSIGNEO](https://img.shields.io/badge/-INSIGNEO-red.svg?colorA=ffffff&colorB=cf2020&logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAAA4AAAAOCAMAAAAolt3jAAABC1BMVEUAAAD%2FAAC%2FAADMGhrRFxfSHh7VHBzJGxvTISHKICDMHx%2FOHR3RJCTTIyPOISHPICDMHR3NIyPNICDRISHQHx%2FQHx%2FOHx%2FQISHOICDQICDQHx%2FOHx%2FPISHQISHQICDOICDQISHOICDQICDPISHPICDQISHPICDPICDOHx%2FPISHPISHQISHPICDPISHPISHQICDOICDOICDPHx%2FPICDPICDOICDPICDPICDPICDPHx%2FPICDPICDQICDPICDPHx%2FQISHOICDPICDPICDPICDPICDPHx%2FMISHOICDPICDQICDPICDPISHPHx%2FPICDPICDPICDPICDPICDPICDOICDPHx%2FPICDPICDPICDPICDi8V76AAAAWHRSTlMAAgQKCxESExcYGRocHR8gIyRITVFSU1ZYYWJjZGZnaWxucXV5fH%2BAg4SFjI%2BUlZeYqKuusbK2ubrDxcbHycvMzdHU2Nrb4uLl5%2Bnq6%2Bzt7u%2F09vf7%2FP3%2B%2FHERCQAAAKVJREFUeAEdx%2BVCg2AAhtHHAANDMQRD7ECwO7DDjViM7b3%2FK9n4zr9Dxdk6jjc4s02si46kK7tMLGDuW%2Fq%2FvdxZlt5shh7VWhsGTqWvTVYll0q4PjMyy7VuYHQ7nIL98oVMAey183sI9EtNPpzXfz7BV8aHIpi%2Fe3ch1iuHysYxnFxHTDf1bD7xpMKBxa6Kk6WVKFXPY8BLZTQWMMYOHvK%2FZHcS6AOapR0V%2FpSSVQAAAABJRU5ErkJggg%3D%3D)](https://insigneo.org/)
[![CompBioMed](https://img.shields.io/badge/-CompBioMed-yellow.svg?colorA=grey&colorB=f4b540&logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAAAoAAAAQCAMAAAAYoR5yAAAA81BMVEUAAAD%2F%2F%2F%2Bqqqr%2F%2F8z%2FzLP%2FyKT%2F1ar%2F4KP%2F9dj%2F5L%2Fm1bP%2F3bv%2Fx4Dw4eHbzJn%2F26%2F52Zn%2FzobQypT%2F6LveyLH%2FyID63d3%2F05P%2Fx3nt19fZzbT%2Fv2Dx0sHk0NDPtX7%2Fxln207D53abEr3f8vFjx1Jz%2FxWP%2FvmL%2FwWT%2Fw0f8vk7t0aj12aLBsnzz0qbTvZritmTAq3vWtF6%2FrIHHr4rJsor6vGPYu3vCsIft17rewJLCq3nw267KtYn%2Fv1vKtJL%2FvFPzz5v547H536P9uE%2FLrG7LsnHGsHbyvmXexIzWq130qkP3y3%2F8tUr3yYD6s0T8sz78tETsoiaWAAAAUXRSTlMAAgMFCg4YGRocHh4gIiMjKCorLS4uNDQ3OT1ISkxPUFFTVldYXV5eYWJkZWdnaGlqamttbW9wcXN1dnh5e35%2Bf3%2BAgYSFio6Sl5ydn6GkqK2fR9KlAAAAd0lEQVQI1zXGRQKCABRAwWdhY3d3t9iJ3d7%2FNC74zmo4H%2FH4AGBxoWIzetKUnDFu321c2vtUrdLx0y5j8oj9u3sPzFL9Ps9K86%2FmSho5FK77sgmglGKamLUA6gH6rkzHDzTcBDfe6BJYO6A4DKdDFmpOQG2PuskfKZ4MqTH%2F64gAAAAASUVORK5CYII%3D)](http://www.compbiomed.eu/)

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) [![Mentioned in Awesome Julia.jl](https://awesome.re/mentioned-badge.svg)](https://github.com/svaksha/Julia.jl/blob/master/Biology.md#bioinformatics)

[![DOI](https://img.shields.io/badge/DOI-10.15131/shef.data.7166183-blue.svg)](https://figshare.com/articles/openBF_Julia_software_for_1D_blood_flow_modelling/7166183)

<a href="https://github.com/INSIGNEO/openBF/actions"><img alt="openbf ci status" src="https://github.com/INSIGNEO/openBF/actions/workflows/ci.yml/badge.svg"></a>

openBF is an open-source 1D blood flow solver based on MUSCL finite-volume numerical scheme, written in [Julia](https://julialang.org/downloads/) and released under [Apache 2.0](http://www.apache.org/licenses/LICENSE-2.0) free software license.

[_Docs_](#docs) | [_Installation_](#installation) | [_Example_](#example) | [_Plot_](#plotting) | [_Dev_](#how-to-dev) | [_Hub_](#ecosystem) | [_Cite_](#cite)

### Docs

- openBF __≦v0.6.3__ check this [website](https://INSIGNEO.github.io/openBF/Docs/index.html) for documentation and tutorials.
- openBF __v0.7+__ docs can be built as `julia make.jl` inside the `docs/` folder (this requires [Documenter](https://github.com/JuliaDocs/Documenter.jl) package).

### Installation

Provided you already have a Julia v1.x installation (all platforms [download](https://julialang.org/downloads/) and Windows [instructions](http://wallyxie.com/weblog/adding-julia-windows-path-command-prompt/)), you can add openBF as

```
julia> ]
(v1.x) pkg> add https://github.com/INSIGNEO/openBF
```

update it as

```
julia> ]
(v1.x) pkg> update openBF
```

test it as

```
julia> ]
(v1.x) pkg> test openBF
```

and use it as

```julia
julia> using openBF
julia> openBF.runSimulation("<input file name>.yml")
```

You can also create (MacOSX/Linux only) an openBF alias as

```bash
$ echo "alias openBF='cp ~/.julia/v1.x/openBF/main.jl ./main.jl && julia main.jl $1'" >> ~/.bashrc
$ source ~/.bashrc
$ openBF -h
usage: main.jl [-v] [-f] [-c] [-h] input_filename

positional arguments:
  input_filename   .yml input file name

optional arguments:
  -v, --verbose    Print STDOUT - default false
  -f, --out_files  Save complete results story rather than only the
                   last cardiac cycle
  -c, --conv_ceil  Ceil convergence value to 100 mmHg (default true)
  -h, --help       show this help message and exit
```

### Example

```
project name: <project name>
results folder: <path/to/your/results/directory>

blood:
  rho: 1060.0 # density
  mu: 4.e-3   # kinematic viscosity

solver:
  Ccfl: 0.9   # Courant number
  cycles: 100 # max number of cardiac cycles
  jump: 100
  convergence tolerance: 15.0 # convergence min error threshold

network:
  - label: <vessel name>
	sn: 1
    tn: 2
    E: 400000 # Young's modulus
    L: 0.04   # length
    R0: 0.012 # lumen radius
    inlet number: 1
    inlet: Q
    inlet file: <inlet filename>
  - label: <vessel name>
  	sn: 2
    tn: 3
    E: 400000
    L: 0.103
    R0: 0.010
    outlet: wk2
  	Cc: 8.2e-11      # peripheral compliance
    R1: 8480000000.0 # peripheral resistance
```

### Plotting

Install [Plots.jl](https://github.com/JuliaPlots/Plots.jl)

```
julia> ]
(v1.x) pkg> add Plots
(v1.x) pkg> <backspace>
julia> using Plots
```
(the first time you run this, the library will be compiled and it may takes several minutes)

plot something
```
using DelimitedFiles

# open the result files
v1 = DelimitedFiles.readdlm("v1_P.last")
v2 = DelimitedFiles.readdlm("v2_P.last")

plot(v1[:,1], v1[:,end]/133.332, label="v1")
plot(v2[:,1], v2[:,end]/133.332, label="v2")

xlabel!("time (s)")
ylabel!("pressure (mmHg)")
```
![waveforms](https://user-images.githubusercontent.com/4661737/97332078-f4101d80-1871-11eb-970d-b7761c069688.png)

### How to dev

```
$ git clone https://github.com/INSIGNEO/openBF.git
$ cd openBF
$ julia
julia> ]
(v1.x) pkg> add Revise
(v1.x) pkg> activate .
(openBF) <backspace>
julia> using Revise
julia> using openBF
```

### Ecosystem

- A collection of 1D networks (from literature) solved by means of openBF can be found in the [openBF-hub](https://github.com/alemelis/openBF-hub) repository.
- The scripts to generate a virtual population of ADAN56s can be found in [openBF-db](https://github.com/alemelis/openBF-db) repository.
- The scripts to generate aged versions of a template openBF model can be found in [openBF_ageing](https://github.com/ibenemerito88/openBF_ageing).
- openBF badge [![openBF](https://img.shields.io/badge/-openBF-red.svg?colorA=ffffff&colorB=008080&logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAABQAAAAOCAQAAACFzfR7AAAA10lEQVQoz4XQIUvDARCG8R%2BCaUEEqwwMwyKCYbImLCp%2BAoPJpMWyYFnfN1BBP4IGk0FkYBkyDYJgkBXTUOcYbFM8g0M3t%2F98rhz3Pnfh6GfKg7bQ8exG1hh2xE%2B9SCeLuT4xlC39RkeKUr1%2B1uWAGMKF%2Be9wW7gzgzX1IS2EhhVIaQnnSj6HlEerCsKrObgeeSeEPfvy7oUzqCaKVzLy3oUnpnWFsi3txIUP6%2BwKdQvIuh2pvdmAUwcyvfdM2PwjNx0mvz3tRAgVOZPGciw0LfqXmprlwdEX%2F8%2BRhjBYrRoAAAAASUVORK5CYII%3D)](https://github.com/INSIGNEO/openBF)

### Publications

openBF has been used in the following works:

- Benemerito I, Narata AP, Narracott A, Marzo A. [Determining clinically-viable biomarkers for ischaemic stroke through a mechanistic and machine learning approach](https://link.springer.com/article/10.1007/s10439-022-02956-7). _Annals of Biomedical Engineering_, 2022. DOI: 10.1007/s10439-022-02956-7
- Mustafa A. [An efficient computational approach to guide intervention in treatment of stroke](https://etheses.whiterose.ac.uk/29992/). _PhD Thesis_, 2021
- Benemerito I, Jordan B, Mustafa A, Marzo A. [Quantification of the effects of ageing, hypertension and atherosclerosis on flow
reversal during a mechanical thrombectomy procedure](https://www.sheffield.ac.uk/insigneo/overview/events/biomedeng-2021-conference), _BioMedEng21_ September 2021, Sheffield, UK _In proceedings_
- Benemerito I, Narata AP, Narracott A, Marzo A. [Pulsatility indices can inform on distal perfusion following ischaemic stroke](https://www.cmbbe-symposium.com/2021/wp-content/uploads/sites/2/2021/09/Program-CMBBE21-A4.qxp_Detailed.pdf), _CMMBE_ September 2021, Online event, _In proceedings_
- Benemerito I, Narata AP, Narracott A, Marzo A. [Identification of biomarkers for perfusion following an ischaemic event](https://cbmc21.vfairs.com/), _CMBE_ September 2021, Online event, _In proceedings_
- Melis A. [Gaussian process emulators for 1D vascular models](http://etheses.whiterose.ac.uk/19175/). _PhD Thesis_, 2017.
- Melis A, Clayton RH, Marzo A. [Bayesian sensitivity analysis of a 1D vascular model with Gaussian process emulators](http://rdcu.be/AqLm). _International Journal for Numerical Methods in Biomedical Engineering_, 2017. DOI: 10.1002/cnm.2882
- Melis A, Clayton RH, Marzo A. [A more efficient approach to perform sensitivity analyses in 0D/1D cardiovascular models](http://www.compbiomed.net/2015/cmbe-proceedings.htm), _CMBE_ July 2015, Cachan, FR. _In proceedings_

Have you used openBF for your research? Let us know!

### Citation

```
@misc{openBF.jl-2018,
title={openBF: Julia software for 1D blood flow modelling}, 
url={https://figshare.com/articles/openBF_Julia_software_for_1D_blood_flow_modelling/7166183/1}, 
DOI={10.15131/shef.data.7166183}, 
abstractNote={
openBF is an open-source 1D blood flow solver based on MUSCL finite-volume numerical scheme, written in Julia and released under Apache 2.0 free software license.

See https://github.com/INSIGNEO/openBF for the git repository and https://insigneo.github.io/openBF/ for the documentation.
}, 
publisher={figshare}, 
author={Melis, Alessandro}, 
year={2018}, 
month={Oct}}
```
