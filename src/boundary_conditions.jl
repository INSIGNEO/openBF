function set_inlet_bc(t::Float64, dt::Float64, v::Vessel, h::Heart)
    v.Q[1] = inlet_from_data(t, h)
    inlet_compatibility!(dt, v)
end

function inlet_from_data(t::Float64, h::Heart)

    # --------------------------------------------------------------------------
    # Functioning
    # --------------------------------------------------------------------------
    # Inlet flow waveform is stored in `h` data structure inside `input_data`
    # matrix. `input_data` has two columns and as many rows as time steps. The
    # first column (`input_data[:,1]`) contains the time variable. The second
    # column (`input_data[:,2]`) contains the flow variable. The two arrays are
    # copied into two variables, `idt` and `idq` respectively.
    # --------------------------------------------------------------------------
    idt = h.input_data[:, 1]
    idq = h.input_data[:, 2]
    # `t_hat` variable tells us how many cardiac cycles have already been
    # simulated. This is calculated with `div(a, b)` which returns `a/b`,
    # truncating to an integer.
    t_hat = div(t, h.cardiac_period)
    # The current time is then updated by subtracting as many cardiac cycles as
    # those already passed. This is done in order to have a reference within the
    # cardiac cycle, i.e. in order to know at what point (*when*) in the cardiac
    # cycle the simulation is running.
    t -= t_hat * h.cardiac_period
    # The waveform is given in a discrete form and usually with a small amount
    # of points along the cardiac cycle. Usually, the numerical $\Delta t$ is
    # small enough to fall between two time steps reported in `idt`. This
    # loop skims all the entries in `idt` and finds the time interval containing
    # the current time `t`. The variable `idx` is initialised to zero for
    # simplicity.
    #
    # <div style="text-align:center">
    # <img src="images/idt_array.png" width="200"></div>
    idx = 0
    for i = 1:length(idt)
        if ((t >= idt[i]) && (t <= idt[i+1]))
            idx = i
            break
        end
    end
    # The flow `qu` is then computed at the current time by a linear
    # interpolation between the two known flows (`idq[idx]` and `idq[idx+1]`)
    # at the ends of the time interval.
    qu = idq[idx] + (t - idt[idx]) * (idq[idx+1] - idq[idx]) / (idt[idx+1] - idt[idx])
    # --------------------------------------------------------------------------
    # Returns:
    # ---------- ---------------------------------------------------------------
    # `qu`       `::Float` inlet flow rate at current time.
    # --------------------------------------------------------------------------
    return qu
end

function riemann_invariants(i::Int64, v::Vessel)
    W1 = v.u[i] - 4 * v.c[i]
    W2 = v.u[i] + 4 * v.c[i]
    W1, W2
end

function inv_riemann_invariants(W1::Float64, W2::Float64)
    u = 0.5 * (W1 + W2)
    c = (W2 - W1) * 0.125
    u, c
end

function inlet_compatibility!(dt::Float64, v::Vessel)

    W11, W21 = riemann_invariants(1, v)
    W12, W22 = riemann_invariants(2, v)

    W11 += (W12 - W11) * (v.c[1] - v.u[1]) * dt / v.dx
    W21 = 2 * v.Q[1] / v.A[1] - W11

    v.u[1], v.c[1] = inv_riemann_invariants(W11, W21)

    v.A[1] = v.Q[1] / v.u[1]
    v.P[1] = pressure(v.A[1], v.A0[1], v.beta[1], v.Pext)

end


function set_outlet_bc(dt::Float64, v::Vessel)
    # A proximal resistance `R1` set to zero means that a reflection coefficient
    if v.R1 == 0.0
        v.P[end] = 2 * v.P[end-1] - v.P[end-2]
        outlet_compatibility!(dt, v)
    else
        wk3!(dt, v)
    end
end

function outlet_compatibility!(dt::Float64, v::Vessel)

    W1M1, W2M1 = riemann_invariants(v.M - 1, v)
    W1M, W2M = riemann_invariants(v.M, v)

    # $$
    #   W_{2M}^{t+1} = W_2^t + \frac{W_{2 M-1}^t - W_{2M}^t}{\Delta x} (u_M^t
    #                   + c_M^t) \Delta t,
    # $$
    # $$
    #   W_{1M}^{t+1} = W_{1M}^0 - R_t \left(W_{2L}^{t+1} - W_{2L}^0 \right).
    # $$
    W2M += (W2M1 - W2M) * (v.u[end] + v.c[end]) * dt / v.dx
    W1M = v.W1M0 - v.Rt * (W2M - v.W2M0)
    # Here the remaining quatities are computed by means of their definitions.
    v.u[end], v.c[end] = inv_riemann_invariants(W1M, W2M)
    v.Q[end] = v.A[end] * v.u[end]

end

# The three element windkessel simulates the perfusion of downstream vessels.
# This 0D model is coupled by [`wk3!`](boundary_conditions.html#wk3!)
# function to 1D model terminal branches via the solution of a Riemann
# problem at the 0D/1D interface.
#
# <div style="text-align:center">
# <img src="images/wk3!.pdf.png" width="250"></div>
#
# $R_1$ is the proximal resistance, $R_2$ is the peripeheral reisistance,
# $C_c$ is the peripheral compliance, $P_c$ is the pressure across the
# peripheral compliance, $P_out$ is the pressure at the artery-vein interface,
# $P_e$ is the pressure at the 0D/1D interface.

# *function* __`wk3!`__
#
# ----------------------------------------------------------------------------
# Parameters:
# ------------- --------------------------------------------------------------
# `dt`          `::Float` $\Delta t$, current time step.
#
# `v`           `::Vessel` vessel data structure.
# ----------------------------------------------------------------------------
# <a name="wk3!"></a>
function wk3!(dt::Float64, v::Vessel)

    # --------------------------------------------------------------------------
    # Functioning [^3]
    # --------------------------------------------------------------------------
    # At capillary level the pressure is assumed to be zero, i.e. $P_{out}=0.$
    # --------------------------------------------------------------------------
    Pout = 0.0
    # The coupling is performed by assuming that an intermediate state
    # $(A^*,u^*)$ generates from $(A_l,u_l)$ (1D outlet) and $(A_r,u_r)$ (0D
    # inlet). This intermediate state must satisfy the windkessel equation
    # $$
    #   A^*u^* \left(1 + \frac{R_1}{R_2}\right) + C_c R_1 \frac{\partial
    #   (A^*u^*)}{\partial t} = \frac{P_e - P_{out}}{R_2} + C_c \frac{\partial
    #   P_e}{\partial t}.
    # $$
    Al = v.A[end]
    ul = v.u[end]
    # $P_c$ is computed at each time step from
    # $$
    #   C_c \frac{\partial P_c}{\partial t} = A^*u^* - \frac{P_c -
    #   P_{out}}{R_2},
    # $$
    # which is discretised numerically with a first-order scheme. $P_c$ is
    # is [initialised](initialise.html#initialiseVessel) to zero.
    v.Pc += dt / v.Cc * (Al * ul - (v.Pc - Pout) / v.R2)
    # @printf "%8.3e\n" v.Pc/v.R2*1e6
    # We consider $\beta$ and $A_0$ to be the same on both sides of the 0D/1D
    # interface. It yields the nonlinear equation
    # $$
    #   \mathcal{F}(A^*) = A^*R_1\left(u_l+4c_l\right)
    #   -4A^*R_1c^* - \frac{\beta}{A_0}\left(
    #   \sqrt{A^*}-\sqrt{A_0}\right) + P_c = 0,
    # $$
    # where $c_l$ and $c^*$ are the wave speeds calculated with $A_l$ and $A^*$,
    # respectively.
    As = Al

    # fun(As) = v.R1*(ul+4*sqrt(3*v.gamma[end]*sqrt(Al)*0.5))*As -
    # 	4*v.R1*sqrt(3*v.gamma[end]*sqrt(As)*0.5)*As -
    # 	(v.Pext + v.beta[end]*(sqrt(As/v.A0[end]) - 1)) + v.Pc

    # dfun(As) = v.R1*( ul +4*sqrt(1.5*v.gamma[end])*(sqrt(sqrt(Al)) - 1.25*sqrt(sqrt(As))) ) -
    # 	v.beta[end]*0.5/(sqrt(As*v.A0[end]))

    ssAl = sqrt(sqrt(Al))
    sgamma = 2 * sqrt(6 * v.gamma[end])
    sA0 = sqrt(v.A0[end])
    bA0 = v.beta[end] / sA0

    fun(As) =
        As * v.R1 * (ul + sgamma * (ssAl - sqrt(sqrt(As)))) -
        (v.Pext + bA0 * (sqrt(As) - sA0)) + v.Pc

    dfun(As) = v.R1 * (ul + sgamma * (ssAl - 1.25 * sqrt(sqrt(As)))) - bA0 * 0.5 / sqrt(As)

    # $\mathcal{F}(A^*)=0$ is solved for $A^*$ with the Newton's method
    # implemented in [__Roots__](https://github.com/JuliaLang/Roots.jl) library.
    # $A^*$ is initialised equal to $A_l$.
    # As = newton(fun, dfun, As)
    As = newtone(fun, dfun, As)
    # Once $A^*$ is found, $u^*$ reads
    # $$
    #   u^* = \frac{P_e^* - P{out}}{A^*R_1},
    # $$
    # where $P_e^*$ is $P_e$ calculated with $A^*$.
    us = (pressure(As, v.A0[end], v.beta[end], v.Pext) - Pout) / (As * v.R1)

    v.A[end] = As
    v.u[end] = us

end

function newtone(f::Function, df::Function, x0)
    xn = x0 - f(x0) / df(x0)
    if abs(xn - x0) <= 1e-5
        return xn
    else
        newtone(f, df, xn)
    end
end

update_ghost_cells!(n::Network) = update_ghost_cells!.(values(n.vessels))
function update_ghost_cells!(v::Vessel)
    v.U00A = v.A[1]
    v.U00Q = v.Q[1]
    v.U01A = v.A[2]
    v.U01Q = v.Q[2]

    v.UM1A = v.A[v.M]
    v.UM1Q = v.Q[v.M]
    v.UM2A = v.A[v.M-1]
    v.UM2Q = v.Q[v.M-1]
end
