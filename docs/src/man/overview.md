# Overview

The cardiovascular system is a network of elastic vessels through which blood is pumped by the heart. The pulsatile flow and vessel elasticity cause pressure to propagate as waves. Bifurcations, bends, and pathological features cause partial wave reflection.

## Why 1D?

The 1D model is derived from the Navier-Stokes equations under three assumptions:

- vessels are long, narrow, elastic tubes;
- blood is incompressible (constant density);
- blood behaves as a Newtonian fluid (constant viscosity).

The resulting system of PDEs captures pulse wave transmission and reflection at a fraction of the computational cost of 3D simulations, making whole-network simulations tractable.

## Numerical scheme

openBF solves the 1D blood flow equations using a second-order MUSCL finite-volume scheme with a predictor-corrector time integration. Junctions (confluences, bifurcations, anastomoses) are resolved by a Newton solver enforcing conservation of mass and total pressure. Outlet boundary conditions are either a reflection coefficient or a Windkessel (two- or three-element) model.

## Acknowledgements

openBF was developed at the [INSIGNEO Institute for in silico Medicine](https://www.sheffield.ac.uk/insigneo), University of Sheffield, as part of A. Melis' PhD thesis under the supervision of [Dr. A. Marzo](https://www.sheffield.ac.uk/mecheng/academic-staff/alberto-marzo). Funding from the UK EPSRC (Grant EP/K037145/1) is gratefully acknowledged.
