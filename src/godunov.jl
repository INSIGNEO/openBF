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

# 1D blood flow within elastic vessels is ruled by a system of hyperbolic
# partial derivative equations
# $$
#   \begin{cases}
#     \dfrac{\partial A}{\partial t} + \dfrac{\partial Au}{\partial x} = 0, \\
#     \dfrac{\partial Au}{\partial t} + \dfrac{\partial Au^2}{\partial x}
#       \dfrac{A}{\rho}\dfrac{\partial P}{\partial x} = - Ru, \\
#     P(x, t) = P_e(x, t) + K_a \left[\left(\dfrac{A}{A_0} \right)^{1/2}
#     -1 \right],
#   \end{cases}
# $$
# By defining $\mathbf{U} = [A, u]^T$, the conservative form of the system
# can be written along with initial and boundary conditions as
# $$
#   \begin{cases}
#     \mathbf{U}_t + \mathbf{F}(\mathbf{U})_x = \mathbf{S}(\mathbf{U}),
#       \quad x \in \big[0, L \big], \\
#     \mathbf{U}(x, 0) = \mathbf{U}^{(0)}(x), \quad t>0, \\
#     \mathbf{U}(0, t) = \mathbf{U}_l(t), \mathbf{U}(L, t) = \mathbf{U}_r(t),
#   \end{cases}
# $$
# where $\mathbf{F}$ is the *flux* term, and $\mathbf{S}$ is the *source*
# term. This set of equations and initial and
# boundary conditions is also called Initial Boundary Value Problem (IBVP).
# The IBVP is solved by means of the Godunov's method [^1], a
# conservative first order finite volume method. The method is applied on
# a computational domain (red) subdivided in $M$ *cell*s $I_i$ (gray).

# <div style="text-align:center">
# <img src="images/comp_domain.pdf.png" width="400"></div>

# Domain length is
# equal to vessel length $L$, and cell length is given by $\Delta x$ size.
# The total simulation time $T$ defines domain height, hence cell height is
# given by $\Delta t$. In each cell $I_i$ three points are defined: center point $x_i$, left
# extreme $x_{i-1/2}$, and right extreme $x_{i+1/2}$.

# <div style="text-align:center">
# <img src="images/cell_i.pdf.png" width="350"></div>

# Godunov's method assumes a piece-wise constant distribution of $\mathbf{U}$
# at each $t = t^n$.

# <div style="text-align:center">
# <img src="images/piecewise.pdf.png" width="350"></div>

# For a *small enough* $\Delta t$, the method updates the solution as
# $$
#   \mathbf{U}_i^{n+1} = \mathbf{U}_i^{n} + \frac{\Delta t}{\Delta x}
#     \left(\mathbf{F}_{i-1/2} - \mathbf{F}_{i+1/2} \right),
# $$
# where $\mathbf{F}_{i\pm 1/2} = \mathbf{F}\left(\mathbf{U}_{i\pm 1/2}(0)
# \right)$ is the *numerical* flux. The argument of the numerical flux,
# $\mathbf{U}_{i\pm 1/2}(0)$, is computed as the solution of a homogeneous (
# the source term will be added later on) quasi-linear local Riemann
# problem (RP) at the interface of two cells. For $\mathbf{U}_{i + 1/2}(0)$
# we have
# $$
#   \begin{cases}
#     \mathbf{U}_t + \mathbf{F}_\mathbf{U} \mathbf{U}_x = 0, \\
#     \mathbf{U}(x_i, 0) = \mathbf{U}^{(0)}(x_i) = \begin{cases}
#       \mathbf{U}_L = \mathbf{U}(x_{i-1/2}), \quad x<0, \\
#       \mathbf{U}_R = \mathbf{U}(x_{i+1/2}), \quad x>0. \\
#     \end{cases}
#   \end{cases}
# $$
# The solution of the RP is depicted below. <a name="RP_solution"></a>

# <div style="text-align:center">
# <img src="images/rp_solution.pdf.png" width="400"></div>

# The solution structure includes two families of waves associated to the
# two system Riemann invariants $\lambda_1$ and $\lambda_2$, respectively.
# There are three constant states $\mathbf{U}_L$, $\mathbf{U}^*$, and
# $\mathbf{U}_R$, of which $\mathbf{U}^*$ is unknown and refers
# to the solution in the star region (yellow). We need *jump*
# conditions to cross the waves and connect these states. Depending on the
# nature of each wave we have *shock* or *rarefaction* conditions.

# [`Godunov`](godunov.html#Godunov) function glues together all the functions
# needed to [solve the RP](godunov.html#riemannSolver)
# at each interface and to compute the [numerical flux](godunov.html#flux).

# *function* __`Godunov`__
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `i`              `::Int` cell index within the vessel numerical domain.
#
# `v`              `::Vessel` current vessel data structure.
#
# `dt`             `::Float` $\Delta t$, current time step.
#
# `b`              `::Blood` blood data structure.
# ----------------------------------------------------------------------------
#
# --------------------------------------------------------------------------
# Functioning:
# --------------------------------------------------------------------------
# In order to solve the local RP, local boundary conditions should be
# retrieved at both left ($i-1/2$) and right ($i+1/2$) interfaces.
# --------------------------------------------------------------------------
# <a name="Godunov"></a>
function Godunov(i :: Int64, v :: Vessel, dt :: Float64, b :: Blood)

  # <div style="text-align:center">
  # <img src="images/ghost_cells.pdf.png" width="400"></div>

  # First we compute the states $\mathbf{U}_L$ and $\mathbf{U}_R$ at the
  # left hand side interface ($i-1/2$). When we deal with the inlet node
  # ($i=1$), blue ghost cells are used to find left states. Right states are
  # always inside the numerical domain and they do not need ghost cells.
  #i-1/2 local left
  if i-1 == 0      #1st ghost cell
    UlA = v.U00A
    UlQ = v.U00Q

  elseif i-1 == -1 #2nd ghost cell
    UlA = v.U01A
    UlQ = v.U01Q

  else
    UlA = v.A[i-1]
    UlQ = v.Q[i-1]
  end

  #local right
  UrA = v.A[i]
  UrQ = v.Q[i]
  # Once the local states have been computed, the solution of the RP at the
  # left
  # interface is found by [`riemannSolver`](godunov.html#riemannSolver)
  # function.
  Al, Ql = riemannSolver(UlA, UlQ, UrA, UrQ, v)

  # The same process is applied to the local right interface.
  #i+1/2 local right
  if i+1 == v.M+1     #1st ghost cell
    UrA = v.UM1A
    UrQ = v.UM1Q

  elseif i+1 == v.M+2 #2nd ghost cell
    UrA = v.UM2A
    UrQ = v.UM2Q

  else
    UrA = v.A[i+1]
    UrQ = v.Q[i+1]
  end

  #local left
  UlA = v.A[i]
  UlQ = v.Q[i]

  Ar, Qr = riemannSolver(UlA, UlQ, UrA, UrQ, v)

  # When both interfaces are solved, left and right ($\mathbf{F}_{i\mp 1/2}$)
  # numerical fluxes are computed by [`flux`](godunov.html#flux) fucntion.
  Fl1, Fl2 = flux(Al, Ql, v.gamma)
  Fr1, Fr2 = flux(Ar, Qr, v.gamma)

  # Conservative variables `A` and `Q` are updated by Godunov's method
  # formula.
  v.A[i] += dt/v.dx * (Fl1-Fr1)
  v.Q[i] += dt/v.dx * (Fl2-Fr2)

  # The source term is solved by [`source`](godunov.html#source) function
  # and `Q` is directly updated. The source term is always zero for the
  # continuity equation, hence `A` is not updated.
  v.Q[i] += source(dt, v.A[i], v.u[i], b)

  # Eventually, primitive variables `P`, `u`, and `c` are updated accordingly
  # to their definitions by [`pressure`](converter.html#pressure) and
  # [`waveSpeed`](converter.html#waveSpeed) functions.
  v.P[i] = pressure(v.A[i], v.A0, v.beta, v.Pext)
  v.u[i] = v.Q[i]/v.A[i]
  v.c[i] = waveSpeed(v.A[i], v.gamma)

end

# The [RP solution](godunov.html#RP_solution) is herein coded. The Riemann
# solver computes the cross-sectional area in the star region
# $A^*$ starting from a guess solution (two rarefactions solution). The
# value of $A^*$ comes from the solution of the nonlinear equation
# $$
#   f(A^*) = f_L(A^*) + f_R(A^*) + \Delta u = 0,
# $$
# where $\Delta u = u_R - u_L$, and $f_L$ and $f_R$ are the local fluxes
# computed by [`fl`](godunov.html#fl) and [`fr`](godunov.html#fr),
# respectively.

# *function* __`riemannSolver`__ $\rightarrow$ `As::Float`, `Qs::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `Al`             `::Float` `A` value at the left side of the local
#                  cell interface.
#
# `Ql`             `::Float` `Q` value at the left side of the local
#                  cell interface.
#
# `Ar`             `::Float` `A` value at the right side of the local
#                  cell interface.
#
# `Qr`             `::Float` `Q` value at the right side of the local
#                  cell interface.
#
# `v`              `::Vessel` current vessel data structure.
# ----------------------------------------------------------------------------
# <a name="riemannSolver"></a>
function riemannSolver(Al :: Float64, Ql :: Float64,
                       Ar :: Float64, Qr :: Float64,
                        v :: Vessel)
  # --------------------------------------------------------------------------
  # Functioning:
  # --------------------------------------------------------------------------
  # Primitives are initialised by their definitions, and the two rarefaction
  # solution is applied to compute initial values for $c^*$ and $A^*$.
  # --------------------------------------------------------------------------
  #initialise solution
  ur = Qr/Ar
  ul = Ql/Al
  cl = waveSpeed(Al, v.gamma)
  cr = waveSpeed(Ar, v.gamma)
  du = ur - ul

#     #vacuum
#   if ul + 4*cl <= ur - 4*cr
# #     println("v")
# #     return (Al+Ar)*0.5, (Ql+Qr)*0.5
#     return v.A0, (Ql+Qr)*0.5
#   end

  #two rarefaction solution
  cs = 0.5 * (cl + cr) - (ur - ul)*0.125
  As = 4/9 * (cs^4)/(v.gamma*v.gamma)
  # The nonlinear equation for the RP is defined and solved by Newton's
  # method.
  fun(As) = fL(As, Al, cs, cl, v) + fR(As, Ar, cs, cr, v) + du

  As = newton(fun, As)

  # Eventually, $u^*$ is computed as
  # $$
  #   \frac{1}{2}\left(u_L + u_R\right) + \frac{1}{2}\left[f_R\left(A^*\right)
  #     -f_L\left(A^* \right) \right],
  # $$
  us = 0.5*((ur+ul) + (fR(As, Ar, cs, cr, v) - fL(As, Al, cs, cl, v)))
  # $Q^*$ is calculated from its definition and returned.
  Qs = As * us
  # --------------------------------------------------------------------------
  # Returns:
  # ---------- ---------------------------------------------------------------
  # `As`       `::Float` cross-sectional area in the star region, solution
  #            of the RP at the interface.
  #
  # `Qs`       `::Float` flow rate in the star region, solution
  #            of the RP at the interface.
  # --------------------------------------------------------------------------
  return As, Qs
end

# Local fluxes (left and right) are computed depending on the nature of the
# wave system located at the interface. In particular
# $$
#   f_L = \begin{cases}
#       4\left(c^* - c_L\right), & A^*\leq A_L \quad \text{(rarefaction)}, \\
#       \sqrt{\dfrac{\gamma\left(A^*-A_L\right)\left(A^{* 3/2}-
#         A_L^{3/2}\right)}{A_L A^*}}, & A^*>A_L \quad \text{(shock)},
#     \end{cases}
# $$
# $$
#   f_R = \begin{cases}
#       4\left(c^* - c_R\right), & A^*\leq A_R \quad \text{(rarefaction)}, \\
#       \sqrt{\dfrac{\gamma\left(A^*-A_R\right)\left(A^{* 3/2}-
#         A_R^{3/2}\right)}{A_R A^*}}, & A^*>A_R \quad \text{(shock)},
#     \end{cases}
# $$

# *function* __`fL`__ $\rightarrow$ `left_flux::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# ------------- --------------------------------------------------------------
# `A`           `::Float` current value of $A^*$ at the left interface.
#
# `Al`          `::Float` cross sectional area of the left cell.
#
# `c`           `::Float` current value of $c^*$ at the left interface.
#
# `cl`          `::Float` wave speed in the left cell.
#
# `v`           `::Vessel` current vessel data structure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ------------- --------------------------------------------------------------
# `left_flux`   `::Float` local flux at the left interface.
# ----------------------------------------------------------------------------
# <a name="fL"></a>
function fL(A, Al :: Float64, c :: Float64, cl :: Float64, v :: Vessel)

  if A <= Al
    return 4*(c - cl)

  else
    return sqrt(
                (v.gamma*(A - Al)*(A^1.5 - Al^1.5)) /
                              (Al*A)
               )
  end

end

# *function* __`fR`__ $\rightarrow$ `right_flux::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# ------------- --------------------------------------------------------------
# `A`           `::Float` current value of $A^*$ at the right interface.
#
# `Ar`          `::Float` cross sectional area of the right cell.
#
# `c`           `::Float` current value of $c^*$ at the right interface.
#
# `cr`          `::Float` wave speed in the right cell.
#
# `v`           `::Vessel` current vessel data structure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ------------- --------------------------------------------------------------
# `right_flux`   `::Float` local flux at the right interface.
# ----------------------------------------------------------------------------
# <a name="fR"></a>
function fR(A, Ar :: Float64, c :: Float64, cr :: Float64, v :: Vessel)

  if A <= Ar
    return 4*(c - cr)

  else
    return sqrt(
                (v.gamma*(A - Ar)*(A^1.5 - Ar^1.5)) /
                              (Ar*A)
               )
  end

end

# The flux in conservative form reads
# $$
#   \mathbf{F}(\mathbf{U}) = \left\{ \begin{array}{c}
#       Au \\
#       Au^2 + \gamma A^{3/2}
#     \end{array} \right\},
# $$
# and its numerical form is computed by [`flux`](godunov.html#flux) function.

# *function* __`flux`__ $\rightarrow$ `f1::Float`, `f2::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# ------------- --------------------------------------------------------------
# `A`           `::Float` cross sectional area along the cell.
#
# `Q`           `::Float` flow rate within the cell.
#
# `gamma`       `::Vessel` elastic constant for `c` calculation.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ---------- -----------------------------------------------------------------
# `f1`       `::Float` first component of the numerical flux.
#
# `f2`       `::Float` second component of the numerical flux.
# ----------------------------------------------------------------------------
# <a name="flux"></a>
function flux(A :: Float64, Q :: Float64, gamma :: Float64)

  f1 = Q
  f2 = Q*Q/A + gamma*A^1.5

  return f1, f2
end

# *function* __`source`__ $\rightarrow$ `S::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# ------------ ----------------------------------------------------------------
# `dt`         `::Float` $\Delta t$, current time step.
#
# `A`          `::Float` cross sectional area along the cell.
#
# `u`          `::Vessel` longitudinal velocity along the cell.
#
# `b`          `::Blood` blood data structure.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# The numerical source is updated by means of a first order implicit backward
# Euler method.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Returns:
# ---------- -----------------------------------------------------------------
# `S`        `::Float` source term.
# ----------------------------------------------------------------------------
# <a name="source"></a>
function source(dt :: Float64, A :: Float64, u :: Float64, b :: Blood)

  return -b.Cf*u*dt

end

# *function* __`calculateDeltaT`__ $\rightarrow$ `DT::Float`
#
# ----------------------------------------------------------------------------
# Parameters:
# ---------------- -----------------------------------------------------------
# `vessels`        `::Array{Vessel, 1}`
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# This function calculates the smallest $\Delta t$ allowed by all vessels
# geometries and properties.
# ----------------------------------------------------------------------------
# <a name="calculateDeltaT"></a>
function calculateDeltaT(vessels, dt :: Array{Float64, 1})

  # For each vessel the $\Delta t$ is computed as
  # $$
  #   \Delta t = \frac{\Delta x}{S_{max}} C_{cfl}
  # $$
  # where $S_{max}$ is the maximum between the forward and the backward
  # characteristics, and $C_{cfl}$ is the Courant-Friedrichs-Lewy condition
  # defined by the user.
  i = 1
  for v in vessels

    lambdap = zeros(Float64, v.M)
    for j in 1:v.M
      lambdap[j] = v.u[j] + v.c[j]
    end

    Smax = maximum(abs.(lambdap))

    dt[i] = v.dx*v.Ccfl/Smax

    i+=1
  end
  #
  # --------------------------------------------------------------------------
  # Returns:
  # ----------- --------------------------------------------------------------
  # `DT`        `::Float` $\Delta t$ computed as the smallest between
  #             the smallest $\Delta t$ associated with each vessel and
  #             the $\Delta t$ limiter $0.001$.
  # --------------------------------------------------------------------------
  return minimum(dt)

end

# *function* __`solveModel`__
#
# ----------------------------------------------------------------------------
# Parameters:
# --------------- ------------------------------------------------------------
# `grafo`         `::GenericGraph` graph data structure containing all the
#                 vessels.
#
# `edgess`        `::Array` collection of `edge` objects.
#
# `vessels`       `::Array{Vessel, 1}` collection of vessels in the system.
#
# `heart`         `::Heart` heart data structure.
#
# `blood`         `::Blood` blood data structure.
#
# `dt`            `::Float` $\Delta t$, current time step.
#
# `current_time`  `::Float` current simulation time.
# ----------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------
# Functioning
# ----------------------------------------------------------------------------
# See the [detailed example](grafo.html).
#
# ----------------------------------------------------------------------------
# <a name="solveModel"></a>
function solveModel(vessels, heart :: Heart, edge_list,
                    blood :: Blood, dt :: Float64, current_time :: Float64)

  # for i in 1:length(edgess)
    # edge = edgess[i]


  # edge_index(e::Edge) = edgemap[e]
  # i = 1
  # for e in LightGraphs.edges(grafo)
  #   s = LightGraphs.src(e)
  #   t = LightGraphs.dst(e)
  #
  #   if LightGraphs.indegree(grafo, s) == 0
  #     openBF.setInletBC(current_time, dt, vessels[i], heart)
  #   end
  #
  #   # *Note*: [MUSCL](MUSCL.html#MUSCL) solver is herein used. To use the first-order
  #   # [Godunov](godunov.html#Godunov) method,
  #   # replace `openBF.MUSCL(vessels[i], dt, blood)` with
  #   #
  #   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ {.julia}
  #   #   for j in 1:vessels[i].M
  #   #     openBF.Godunov(j, vessels[i], dt, blood)
  #   #   end
  #   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   openBF.MUSCL(vessels[i], dt, blood)
  #   if LightGraphs.outdegree(grafo, t) == 0
  #     openBF.setOutletBC(dt, vessels[i])
  #
  #   elseif LightGraphs.outdegree(grafo, t) == 1
  #     if length(LightGraphs.in_neighbors(grafo, t)) == 1
  #       o = LightGraphs.out_neighbors(grafo, t)
  #       openBF.joinVessels(blood, vessels[i], vessels[edgemap[o]])
  #
  #     else
  #       es = LightGraphs.in_neighbors(grafo, t)
  #       if i == max(Graphs.edge_index(es[1], grafo), Graphs.edge_index(es[2], grafo))
  #         a = Graphs.edge_index(es[1], grafo)
  #         b = Graphs.edge_index(es[2], grafo)
  #         c = Graphs.edge_index(LightGraphs.out_neighbors(grafo, t)[1])
  #         openBF.solveAnastomosis(vessels[a], vessels[b], vessels[c])
  #       end
  #     end
  #
  #   elseif LightGraphs.outdegree(grafo, t) == 2
  #     openBF.joinVessels(blood, vessels[i], vessels[Graphs.edge_index(Graphs.out_edges(t,grafo)[1])],
  #                       vessels[Graphs.edge_index(Graphs.out_edges(t,grafo)[2])])
  #   end
  # end
  # i += 1
  for j in 1:size(edge_list)[1]
    i = edge_list[j,1]
    s = edge_list[j,2]
    t = edge_list[j,3]
    v = vessels[i]
    lbl = v.label

    if size(find(edge_list[:,3] .== s))[1] == 0
      openBF.setInletBC(current_time, dt, v, heart)
    end

    openBF.MUSCL(v, dt, blood)

    if size(find(edge_list[:,2] .== t))[1] == 0
      openBF.setOutletBC(dt, v)
      # println("\t Outlet vessel - Compute outlet BC")

    elseif size(find(edge_list[:,2] .== t))[1] == 2
      d1_i = find(edge_list[:,2] .== t)[1]
      d2_i = find(edge_list[:,2] .== t)[2]

      openBF.joinVessels(blood, v, vessels[d1_i], vessels[d2_i])
      # println("\t\t Bifurcation at node $t between vessels $i, $d1_i, and $d2_i")

    elseif size(find(edge_list[:,3] .== t))[1] == 1
      d_i = find(edge_list[:,2] .== t)[1]

      openBF.joinVessels(blood, v, vessels[d_i])
      # println("\t\t Junction at node $t between vessels $i and $d_i")

    elseif size(find(edge_list[:,3] .== t))[1] == 2
      p1_i = find(edge_list[:,3] .== t)[1]
      p2_i = find(edge_list[:,3] .== t)[2]
      if maximum([p1_i, p2_i]) == i
        p2_i = minimum([p1_i, p2_i])
        d = find(edge_list[:,2] .== t)[1]
        openBF.joinVessels(blood, v, vessels[p2_i], vessels[d])
        # println("\t\t Anastomosis at node $s between vessels $i, $p2_i, and $d")
      end
    end
  end
end

# ### References
#
# [^1]: Toro, Eleuterio F. Riemann solvers and numerical methods for fluid
# dynamics: a practical introduction. Springer Science & Business Media, 2009.
