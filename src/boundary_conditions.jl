#=
Copyright 2015-2024 INSIGNEO Institute for in silico Medicine

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

function set_inlet_bc(t::Float64, dt::Float64, v::Vessel, h::Heart)
    v.Q[1] = inlet_from_data(t, h)
    inlet_compatibility(dt, v)
end

function inlet_from_data(t::Float64, h::Heart)
    idt = h.input_data[:, 1]
    idq = h.input_data[:, 2]
    t_hat = div(t, h.cardiac_period)
    t -= t_hat * h.cardiac_period
    idx = 0
    for i = 1:length(idt)
        if ((t >= idt[i]) && (t <= idt[i+1]))
            idx = i
            break
        end
    end
    idq[idx] + (t - idt[idx]) * (idq[idx+1] - idq[idx]) / (idt[idx+1] - idt[idx])
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

function inlet_compatibility(dt::Float64, v::Vessel)

    W11, W21 = riemann_invariants(1, v)
    W12, W22 = riemann_invariants(2, v)

    W11 += (W12 - W11) * (v.c[1] - v.u[1]) * dt / v.dx
    W21 = 2 * v.Q[1] / v.A[1] - W11

    v.u[1], v.c[1] = inv_riemann_invariants(W11, W21)

    v.A[1] = v.Q[1] / v.u[1]
    v.P[1] = pressure(v.A[1], v.A0[1], v.beta[1], v.Pext)

end


function set_outlet_bc(dt::Float64, v::Vessel, ρ::Float64)
    # A proximal resistance `R1` set to zero means that a reflection coefficient
    if v.R1 == 0.0
        v.P[end] = 2 * v.P[end-1] - v.P[end-2]
        outlet_compatibility(dt, v)
    else
        wk3(dt, v, ρ)
    end
end

function outlet_compatibility(dt::Float64, v::Vessel)

    W1M1, W2M1 = riemann_invariants(v.M - 1, v)
    W1M, W2M = riemann_invariants(v.M, v)

    W2M += (W2M1 - W2M) * (v.u[end] + v.c[end]) * dt / v.dx
    W1M = v.W1M0 - v.Rt * (W2M - v.W2M0)

    v.u[end], v.c[end] = inv_riemann_invariants(W1M, W2M)
    v.Q[end] = v.A[end] * v.u[end]
end


function wk3(dt::Float64, v::Vessel, ρ::Float64)
    Pout = 0.0
    Al = v.A[end]
    ul = v.u[end]

    # inlet impedance matching
    if v.inlet_impedance_matching
        v.R1 = ρ * wave_speed(v.A[end], v.gamma[end]) / v.A[end]
        v.R2 = v.total_peripheral_resistance - v.R1
    end

    v.Pc += dt / v.Cc * (Al * ul - (v.Pc - Pout) / v.R2)
    As = Al

    ssAl = sqrt(sqrt(Al))
    sgamma = 2 * sqrt(6 * v.gamma[end])
    sA0 = sqrt(v.A0[end])
    bA0 = v.beta[end] / sA0

    fun(As) =
        As * v.R1 * (ul + sgamma * (ssAl - sqrt(sqrt(As)))) -
        (v.Pext + bA0 * (sqrt(As) - sA0)) + v.Pc

    dfun(As) = v.R1 * (ul + sgamma * (ssAl - 1.25 * sqrt(sqrt(As)))) - bA0 * 0.5 / sqrt(As)
    As = newtone(fun, dfun, As)
    us = (pressure(As, v.A0[end], v.beta[end], v.Pext) - Pout) / (As * v.R1)

    v.A[end] = As
    v.u[end] = us
end

# function newtone(f::Function, df::Function, x0)
#     xn = x0 - f(x0) / df(x0)
#     if abs(xn - x0) <= 1e-5
#         return xn
#     else
#         newtone(f, df, xn)
#     end
# end

function newtone(f::Function, df::Function, xn)
    for _=1:10
        xn -= f(xn) / df(xn)
    end
    xn
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
