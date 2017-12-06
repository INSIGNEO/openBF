#= inlet bc selector
3: from data
=#
const inlet_BC_switch = 3
const inlet_type = "Q"

# numerical domain
const Ccfl   = 0.9         # Courant number
const cycles = 100

# blood properties
const rho = 1060.    # density [kg/m3]
const mu  = 4.e-3   # dynamic viscosity [Pa⋅s]
const gamma_profile = 9

# vessel properties
#const Pext = 0.			     # external pressure [Pa]
const initial_pressure = 0.
