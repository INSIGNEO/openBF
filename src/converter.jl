#=
Copyright 2017 INSIGNEO Institute for in silico Medicine

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=#

# This file contains all the functions used within `openBF` to convert
# quantities. In particular, the following functions are provided
#
# *   [`waveSpeed`](converter.html#waveSpeed): $A \rightarrow c$;
# *   [`pressure`](converter.html#pressure): $A \rightarrow P$;
# *   [`areaFromPressure`](converter.html#areaFromPressure):
#       $P \rightarrow A$;
# *   [`riemannInvariants`](converter.html#riemannInvariants):
#       $(u, c) \rightarrow (W_1, W_2)$;
# *   [`rI2uc`](converter.html#rI2uc): $(W_1, W_2) \rightarrow (u, c)$.
#

# Wave speed
# $$
#   c = \sqrt{\frac{3}{2}\gamma \sqrt{A}}, \quad
#     \gamma = \frac{\beta}{3 \rho R_0 \sqrt{\pi}}.
# $$

# *function* __`waveSpeed`__ $\rightarrow$ `c::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `A`              `::Float` cross sectional area.
#
# `gamma`          `::Float` elastic constant.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ------------- --------------------------------------------------------------
# `c`           `::Float` wave speed given cross sectional area A.
# ----------------------------------------------------------------------------
# <a name="waveSpeed"></a>
function waveSpeed(A :: Float64, gamma :: Float64)

  return sqrt(3*gamma*sqrt(A)*0.5)

end

# *function* __`waveSpeed`__ $\rightarrow$ `c::Array{Float, 1}`
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `A`              `::Array{Float, 1}` cross sectional area.
#
# `gamma`          `::Float` elastic constant.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# To compute the wave speed along the entire vessel, `waveSpeed` is
# re-defined to handle also a vector of `A`s. `waveSpeed` is called
# recursively with the same `gamma` along the vessel.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ------------- --------------------------------------------------------------
# `c`           `::Arrary{Float, 1}` wave speed given cross sectional area A.
# ----------------------------------------------------------------------------
# <a name="waveSpeed2"></a>
function waveSpeed(A :: Array{Float64, 1}, gamma :: Array{Float64, 1}, c :: Array{Float64, 1})

  for i in 1:length(A)
    c[i] = waveSpeed(A[i], gamma[i])
  end

  return c
end

# Trans-mural pressure
# $$
#   P = P_{ext} + \beta \left[\left(\frac{A}{A_0}\right)^{1/2} -1 \right], \quad
#     \beta = \sqrt{\frac{\pi}{A_0}} \frac{h_0 E}{1 - \sigma^2}.
# $$

# *function* __`pressure`__ $\rightarrow$ `P::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `A`              `::Float` cross sectional area.
#
# `A0`             `::Float` unstressed cross sectional area.
#
# `beta`           `::Float` elastic constant.
#
# `Pext`           `::Float` constant external pressure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ------------- --------------------------------------------------------------
# `P`           `::Float` trans-mural pressure given the cross sectional area
#               `A`.
# ----------------------------------------------------------------------------
# <a name="pressure"></a>
function pressure(A    :: Float64, A0   :: Float64,
                  beta :: Float64, Pext :: Float64)

  return Pext + beta*(sqrt(A/A0) - 1.)

end

# *function* __`pressure`__ $\rightarrow$ `P::Array{Float, 1}`
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `A`              `::Array{Float, 1}` cross sectional area.
#
# `A0`             `::Float` unstressed cross sectional area.
#
# `beta`           `::Float` elastic constant.
#
# `Pext`           `::Float` constant external pressure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# To compute the pressure along the entire vessel, `pressure` is
# re-defined to handle also a vector of `A`s. `pressure` is called
# recursively with the same `beta` along the vessel.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ------------- --------------------------------------------------------------
# `P`           `::Array{Float, 1}` trans-mural pressure given
#               the cross sectional area `A`.
# ----------------------------------------------------------------------------
# function pressure(A    :: Array{Float64, 1}, A0   :: Float64,
#                   beta :: Float64,           Pext :: Float64)

#   p = zeros(Float64, length(A))

#   for i in 1:length(A)
#     p[i] = Pext + beta*(sqrt(A[i]/A0) - 1.)

#   end

#   return p
# end

function pressure(A    :: Array{Float64, 1}, A0   :: Array{Float64, 1},
                  beta :: Array{Float64,1},  Pext :: Float64,
                  p :: Array{Float64, 1})

  for i in 1:length(A)
    p[i] = Pext + beta[i]*(sqrt(A[i]/A0[i]) - 1.)
  end

  return p
end

# Trans-mural pressure can be inverted to find cross sectional area as
# $$
#   A = A_0 \left( \frac{P-P_{ext}}{\beta} +1 \right)^2.
# $$

# *function* __`areaFromPressure`__ $\rightarrow$ `A::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `P`              `::Float` cross sectional area.
#
# `A0`             `::Float` unstressed cross sectional area.
#
# `beta`           `::Float` elastic constant.
#
# `Pext`           `::Float` constant external pressure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ------------- --------------------------------------------------------------
# `A`           `::Float` cross sectional area given the trans-mural pressure
#               `P`.
# ----------------------------------------------------------------------------
# <a name="areaFromPressure"></a>
function areaFromPressure(P    :: Float64, A0   :: Float64,
                          beta :: Float64, Pext :: Float64)

   return A0 * ((P-Pext)/beta + 1)*((P-Pext)/beta + 1)

end

# *function* __`areaFromPressure`__ $\rightarrow$ `A::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `P`              `::Float` cross sectional area.
#
# `A0`             `::Float` unstressed cross sectional area.
#
# `beta`           `::Float` elastic constant.
#
# `Pext`           `::Float` constant external pressure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# To compute the cross sectional area along the entire vessel,
# `areaFromPressure` is re-defined to handle also a vector of `P`s.
# `areaFromPressure` is called recursively with the same `beta` along the
# vessel.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ------------- --------------------------------------------------------------
# `A`           `::Float` cross sectional area given the trans-mural pressure
#               `P`.
# ----------------------------------------------------------------------------
# function areaFromPressure(P    :: Array{Float64, 1}, A0   :: Float64,
#                           beta :: Float64,           Pext :: Float64)

#   a = zeros(Float64, length(P))

#   for i in 1:length(P)
#     a[i] = areaFromPressure(P[i], A0, beta)

#   end

#   return a
# end

# Riemann invariants are computed from `u` and `c` values as
# $$
#   W_1 = u - 4c, \quad W_2 = u + 4c.
# $$

# *function* __`riemannInvariants`__ $\rightarrow$ `W1::Float`, `w2::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `i`              `::Float` cell index.
#
# `v`              `::Vessel` vessel data structure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ------------- --------------------------------------------------------------
# `W1`          `::Float` backward Riemann invariant.
#
# `W2`          `::Float` forward Riemann invariant.
# ----------------------------------------------------------------------------
# <a name="riemannInvariants"></a>
function riemannInvariants(i :: Int64, v :: Vessel)

  W1 = v.u[i] - 4*v.c[i]
  W2 = v.u[i] + 4*v.c[i]

  return W1, W2
end

# `u` and `c` can be computed from Riemann invariants as
# $$
#   u = \frac{1}{2} (W_1 + w_2), \quad c = \frac{W_2 - W_1}{8}.
# $$

# *function* __`rI2uc`__ $\rightarrow$ `W1::Float`, `w2::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `i`              `::Float` cell index.
#
# `v`              `::Vessel` vessel data structure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ------------- --------------------------------------------------------------
# `u`           `::Float` longitudinal velocity.
#
# `c`           `::Float` wave speed.
# ----------------------------------------------------------------------------
# <a name="rI2uc"></a>
function rI2uc(W1 :: Float64, W2 :: Float64)
  u = 0.5*(W1 + W2)
  c = (W2 - W1)*0.125

  return u, c
end

# *function* __`calculatePrimitives`__
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `v`              `::Vessel` vessel data structure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# This function updates all primitives variables at once.
# ----------------------------------------------------------------------------
# <a name="calculatePrimitives"></a>
# function calculatePrimitives(v :: Vessel)

#   v.P = pressure(v.A, v.A0, v.beta)
#   v.u = v.Q./v.A
#   v.c = waveSpeed(v.A, v.gamma)

# end
