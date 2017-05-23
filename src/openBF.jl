#= Copyright (C) 2017 Alessandro Melis.

  This file is part of openBF.

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Bash.  If not, see <http://www.gnu.org/licenses/>.
=#

# __openBF__ is a computational library written
# in [Julia](julialang.org) for blood flow simulations
# in straight elastic arteries. Navier-Stokes equations are reduced
# to 1D form and solved along the longitudinal direction. Numerical
# solution is achieved via a first-order finite volume solver
# based on the Godunov' method. For an usage example refer to
# [main.jl](main.html) file and its documentation.
module openBF

# ### Import external libraries

# [__Graphs__](https://github.com/JuliaLang/Graphs.jl) "is a
# Julia package that provides graph types and algorithms."
# It is used for arterial tree parsing and described in
# [utilities.jl](godunov.html#solveModel).
using Graphs

# [__ProgressMeter__](https://github.com/timholy/ProgressMeter.jl)
# provides a simple implementation of a loading bar to be shown while
# running the solution.
using ProgressMeter

# ### Import openBF' files

# OpenBF' own types are contained in [BTypes.jl](BTypes.html) where
# data structures for vessels, blood, and numerical scheme
# are specified.
using  BTypes

# Data structures and output files are initialised at the beginning
# of each simulation. [initialise.jl](initialise.html) contains all
# the functions needed to read model description and the
# global settings.
include("initialise.jl")

# Inlet and outlets boundary conditions are applied at each time step,
# and vessels extremity nodes are updated with functions defined in
# [boundary_conditions.jl](boundary_conditions.html).
include("boundary_conditions.jl")

# Interface problems (junctions and bifurcations) are solved via the
# method of characteristics. The function in
# [junctions.jl](junctions.html) discriminates between
# [conjunctions](conjunctions.html) and
# [bifurcations](bifurcations.jl), and calls the appropriate solver.
include("junctions.jl")
include("conjunctions.jl")
include("bifurcations.jl")

# [godunov.jl](godunov.html) contains Godunov' method functions.
include("godunov.jl")

# [MUSCL.jl](MUSCL.html) is where the MUSCL method is implemented.
include("MUSCL.jl")

# Conversions from pressure to cross-sectional area (and vice versa),
# or from cross-sectional area to wave speed (and vice versa) are
# handled by functions in [converter.jl](converter.html); Riemann
# invariants are computed in the same location as well.
include("converter.jl")

# Input and output writing functions are stored defined in
# [IOutils.jl](IOutils.html).
include("IOutils.jl")

# Convergence check functions are all defined in [check_convergence.jl](check_convergence.html).
include("check_convergence.jl")

include("anastomosis.jl")

end
