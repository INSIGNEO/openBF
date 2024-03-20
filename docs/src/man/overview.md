The cardiovascular system is a complex network of elastic vessels through which blood is pumped by contraction of the heart. This pulsatile regime and vessel elasticity cause pressure to propagate along the arterial circulation as waves. Mechanical discontinuities caused by bifurcation, bends or cardiovascular pathology cause pressure waves to be reflected in all directions.

![](https://miro.medium.com/v2/resize:fit:1400/format:webp/1*XFQ5yKLnyEci7hLaQh0nEQ.png)

Pressure waves measured at a specific location can be seen as the result of a superimposition of incident (forward) and reflected (backward) waves. The analysis of this superposition mechanism allows for the study of mechanical properties upstream and downstream of the measurement point and can be a rich source of diagnostic information about the system through which these waves propagate. However, given the vascular system complexity, it is at present difficult to ascribe a particular waveform feature to a specific trait of the arterial circulation

The understanding of the cause-effect mechanisms governing wave propagation in the cardiovascular system was used at [INSIGNEO](https://www.sheffield.ac.uk/insigneo) to develop a computer model that can reproduce this behaviour in a quantitative way. openBF is a Julia package for the solution of one-dimensional (1D) blood flow in networks of elastic arteries.

## Why 1D?

The 1D model is derived from the three-dimensional formulation of the viscous fluid motion equations, the Navier-Stokes equations. These are simplified by assuming the arteries to be long, narrow, and elastic tubes and the blood to be incompressible (i.e., its density is constant over time) and behaving as a Newtonian fluid (i.e., its viscosity does not depend on its velocity).

![](https://miro.medium.com/v2/resize:fit:1400/format:webp/1*z7oY7GOWzD9YsTRiZjxpYw.png)


## Does it work?

Despite all the above assumptions, the resulting partial differential equation system is complex enough to simulate the physics of pulse wave transmission and reflection. This is at the expense of results accuracy in areas where the flow structure is inherently three-dimensional, e.g., in the proximity of heart valves and bifurcations. Conversely, the computational requirements are small compared to 3D simulations, allowing the simulation of the entire cardiovascular system rather than only circumscribed parts.

![](https://miro.medium.com/v2/resize:fit:2000/format:webp/1*HpOB3AYJ2dfPJGiF9Vf6Fw.png)


## Acknowledgements

The software openBF was developed at INSIGNEO Institute for in silico Medicine at The University of Sheffield (UK) as part of A. Melis’ PhD thesis under the supervision of [Dr. A. Marzo](https://www.sheffield.ac.uk/mecheng/academic-staff/alberto-marzo). We gratefully acknowledge funding from the UK Engineering and Physical Sciences Research Council (Grant Number EP/K037145/1).
openBF is based on MUSCL finite-volume numerical scheme, written in Julia, and released under [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0) free software license on GitHub.

The solver is currently used in several HPC systems within the H2020 [CompBioMed](https://www.compbiomed.eu) Centre of Excellence in Computational Biomedicine (Grant Agreement №675451). In [SURFSara](https://www.surf.nl/en/research-it)’s HPC-Cloud, it is being used for large-scale sensitivity analysis and uncertainty quantification studies, towards clinical application. The simulations reported in the thesis and in this [paper](https://pubmed.ncbi.nlm.nih.gov/28337862/) where run in TUoS’s ShARC tier-3 cluster. More information can be found on the official page and in CompBioMed [Software Hub](https://www.compbiomed.eu/services/software-hub/compbiomed-software-openbf/).
