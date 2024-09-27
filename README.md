# __openBF.jl__

[![TUoS](https://img.shields.io/badge/-The%20University%20of%20Sheffield-blue.svg?colorA=ffffff&colorB=009fe3&logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAAA4AAAAOCAMAAAAolt3jAAABsFBMVEUAAABmZplVVYBLPHheUYZJPXlSR3pCOnNgWIdoYI9FPnVMQXc0JGM0KmhEOXI0MW8xMW2Ph6yYkrOWj6%2BOiayinLqWkbIcWJkdWJezrsW2sciclrWmoL0cWpkbXJuCfqQVba4Vba2tqsKBgKWXkbGgm7ioo7%2Bln7umobyppcClob0RdrgSd7ikoLuXkbEeeKwieqmYkrKkn7sAneAAn%2BMBnN8Dn%2BIEj9IEkNQEoOMFi88GicwHnuAInuAVaqkVgrEWa6sWiskYj84ek88mk8wmqOInqOIqhrA4qd45qd87kaA8q988q%2BA9qN0%2Fk51Aqd1CqdxCqd1JnqtNnc9OnKVPnItSmsZUmqZVlsFVnYZVns1bst9ds99ftN9isNtitN9joa1msdtnveZqoaRqrNZqvuZsrdZvk2lylWR0nKZ2oJx%2BfKOBf6WEn5eGr3iIr3eduNWfmbifmrifudWim6mjnKmloLymrz2psTuuqb6uqsOwq8C0taG2tqS4tqy4tra6uTC6ui6%2FsCK%2FsCPBvdHDv9LFwtfGwtjOy9rPzdvUlaHmkpb%2B%2FPz%2B%2Fv7%2F%2FPz%2F%2F%2F9OZpcfAAAAM3RSTlMABQYRExUZHyAgJS8xMTFTVFleZmiJkpOTmp2oqba9wsTFxcfO0dPX3OLj5%2Bf19v39%2Fv4kncL4AAAAvElEQVR4AWNgYBIQV9Q20NdRkOBnZGBg4DExTQ9taIzMMDbhA3PdzMqbmyosXaHc4kCvvHy%2FgAIINzElvL%2B%2FvzcoLQ7E5bX1TOhq624P8bDhBnLZTEwye%2Fr7%2BpNNTFgYgMApJreuprY%2BK94RyOGSlEv19XZx9vFPkpViZ1BRZ7WKdq%2Bsto%2B1ZtZQYhAp1FS2iCopDTaX1yoSZmDgkNHrsAuLcOjUleZkAAGxqtbsnJYyUQYoEFIzMlQVBLEA%2FZgsl9iPrB4AAAAASUVORK5CYII%3D)](https://www.sheffield.ac.uk)
[![INSIGNEO](https://img.shields.io/badge/-INSIGNEO-red.svg?colorA=ffffff&colorB=cf2020&logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAAA4AAAAOCAMAAAAolt3jAAABC1BMVEUAAAD%2FAAC%2FAADMGhrRFxfSHh7VHBzJGxvTISHKICDMHx%2FOHR3RJCTTIyPOISHPICDMHR3NIyPNICDRISHQHx%2FQHx%2FOHx%2FQISHOICDQICDQHx%2FOHx%2FPISHQISHQICDOICDQISHOICDQICDPISHPICDQISHPICDPICDOHx%2FPISHPISHQISHPICDPISHPISHQICDOICDOICDPHx%2FPICDPICDOICDPICDPICDPICDPHx%2FPICDPICDQICDPICDPHx%2FQISHOICDPICDPICDPICDPICDPHx%2FMISHOICDPICDQICDPICDPISHPHx%2FPICDPICDPICDPICDPICDPICDOICDPHx%2FPICDPICDPICDPICDi8V76AAAAWHRSTlMAAgQKCxESExcYGRocHR8gIyRITVFSU1ZYYWJjZGZnaWxucXV5fH%2BAg4SFjI%2BUlZeYqKuusbK2ubrDxcbHycvMzdHU2Nrb4uLl5%2Bnq6%2Bzt7u%2F09vf7%2FP3%2B%2FHERCQAAAKVJREFUeAEdx%2BVCg2AAhtHHAANDMQRD7ECwO7DDjViM7b3%2FK9n4zr9Dxdk6jjc4s02si46kK7tMLGDuW%2Fq%2FvdxZlt5shh7VWhsGTqWvTVYll0q4PjMyy7VuYHQ7nIL98oVMAey183sI9EtNPpzXfz7BV8aHIpi%2Fe3ch1iuHysYxnFxHTDf1bD7xpMKBxa6Kk6WVKFXPY8BLZTQWMMYOHvK%2FZHcS6AOapR0V%2FpSSVQAAAABJRU5ErkJggg%3D%3D)](https://insigneo.org/)
[![CompBioMed](https://img.shields.io/badge/-CompBioMed-yellow.svg?colorA=grey&colorB=f4b540&logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAAAoAAAAQCAMAAAAYoR5yAAAA81BMVEUAAAD%2F%2F%2F%2Bqqqr%2F%2F8z%2FzLP%2FyKT%2F1ar%2F4KP%2F9dj%2F5L%2Fm1bP%2F3bv%2Fx4Dw4eHbzJn%2F26%2F52Zn%2FzobQypT%2F6LveyLH%2FyID63d3%2F05P%2Fx3nt19fZzbT%2Fv2Dx0sHk0NDPtX7%2Fxln207D53abEr3f8vFjx1Jz%2FxWP%2FvmL%2FwWT%2Fw0f8vk7t0aj12aLBsnzz0qbTvZritmTAq3vWtF6%2FrIHHr4rJsor6vGPYu3vCsIft17rewJLCq3nw267KtYn%2Fv1vKtJL%2FvFPzz5v547H536P9uE%2FLrG7LsnHGsHbyvmXexIzWq130qkP3y3%2F8tUr3yYD6s0T8sz78tETsoiaWAAAAUXRSTlMAAgMFCg4YGRocHh4gIiMjKCorLS4uNDQ3OT1ISkxPUFFTVldYXV5eYWJkZWdnaGlqamttbW9wcXN1dnh5e35%2Bf3%2BAgYSFio6Sl5ydn6GkqK2fR9KlAAAAd0lEQVQI1zXGRQKCABRAwWdhY3d3t9iJ3d7%2FNC74zmo4H%2FH4AGBxoWIzetKUnDFu321c2vtUrdLx0y5j8oj9u3sPzFL9Ps9K86%2FmSho5FK77sgmglGKamLUA6gH6rkzHDzTcBDfe6BJYO6A4DKdDFmpOQG2PuskfKZ4MqTH%2F64gAAAAASUVORK5CYII%3D)](http://www.compbiomed.eu/)

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) [![Mentioned in Awesome Julia.jl](https://awesome.re/mentioned-badge.svg)](https://github.com/svaksha/Julia.jl/blob/master/Biology.md#bioinformatics)

[![DOI](https://zenodo.org/badge/92175375.svg)](https://doi.org/10.5281/zenodo.13850604)

<a href="https://github.com/INSIGNEO/openBF/actions"><img alt="openbf ci status" src="https://github.com/INSIGNEO/openBF/actions/workflows/ci.yml/badge.svg"></a>

[![chat](https://dcbadge.vercel.app/api/server/eVU5bDEm8X)](https://discord.gg/eVU5bDEm8X)

[![openBF](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://openBF.streamlit.app)

openBF is an open-source 1D blood flow solver based on MUSCL finite-volume numerical scheme, written in [Julia](https://julialang.org/downloads/) and released under [Apache 2.0](http://www.apache.org/licenses/LICENSE-2.0) free software license.

Read the [documentation](https://insigneo.github.io/openBF/stable) for installation and run instructions.

## Release notes

### v2.0.0
This is a complete re-write of openBF solver with several bugfixes and few new functionalities.
Config files should be backward compatible but please refer to the new [documentation](https://insigneo.github.io/openBF/stable) for more details.

Currently not supported:
- multiple inlets

If your workflow relies on this feature, we recommend to use [release v1.5.1](https://github.com/INSIGNEO/openBF/releases/tag/v1.5.1).

---

### Citation

```
@misc{openBF.jl-2018,
title={openBF: Julia software for 1D blood flow modelling}, 
url={https://figshare.com/articles/openBF_Julia_software_for_1D_blood_flow_modelling/7166183/1}, 
DOI={10.15131/shef.data.7166183}, 
abstractNote={
openBF is an open-source 1D blood flow solver based on MUSCL finite-volume numerical scheme, written in Julia and released under Apache 2.0 free software license.

See https://github.com/INSIGNEO/openBF for the git repository and https://insigneo.github.io/openBF/stable for the documentation.
}, 
publisher={figshare}, 
author={Melis, Alessandro}, 
year={2018}, 
month={Oct}}
```
