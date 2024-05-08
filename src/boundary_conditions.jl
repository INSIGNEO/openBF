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

function inbc!(v::Vessel, t::Float64, dt::Float64, h::Heart)
    v.Q[1] = inlet_from_data(t, h)
    incompat!(v, dt)
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
    c = wave_speed(v.A[i], v.gamma[i+1])
    (v.u[i] - 4c, v.u[i] + 4c)
end

function inv_riemann_invariants(W1::Float64, W2::Float64)
    0.5 * (W1 + W2)
end

function incompat!(v::Vessel, dt::Float64)
    W11, W21 = riemann_invariants(1, v)
    W12, W22 = riemann_invariants(2, v)

    W11 += (W12 - W11) * (wave_speed(v.A[1], v.gamma[2]) - v.u[1]) * dt / v.dx
    W21 = 2 * v.Q[1] / v.A[1] - W11

    v.u[1] = inv_riemann_invariants(W11, W21)
    v.A[1] = v.Q[1] / v.u[1]
end


function outbc!(v::Vessel, dt::Float64, ρ::Float64)
    if v.usewk3
        wk3!(v, dt, ρ)
    else
        outcompat!(v, dt)
    end
end

function outcompat!(v::Vessel, dt::Float64)
    W1M1, W2M1 = riemann_invariants(v.M - 1, v)
    W1M, W2M = riemann_invariants(v.M, v)

    W2M += (W2M1 - W2M) * (v.u[end] + wave_speed(v.A[end], v.gamma[end-1])) * dt / v.dx
    W1M = v.W1M0 - v.Rt * (W2M - v.W2M0)

    v.u[end] = inv_riemann_invariants(W1M, W2M)
    v.Q[end] = v.A[end] * v.u[end]
end


function wk3!(v::Vessel, dt::Float64, ρ::Float64)
    # inlet impedance matching
    if v.inlet_impedance_matching
        v.R1 = ρ * wave_speed(v.A[end], v.gamma[end-1]) / v.A[end]
        v.R2 = abs(v.total_peripheral_resistance - v.R1)
    end

    v.Pc += dt / v.Cc * (v.A[end] * v.u[end] - (v.Pc - v.Pout) / v.R2)
    As = v.A[end]

    ssAl = sqrt(sqrt(v.A[end]))
    sgamma = 2 * sqrt(6 * v.gamma[end-1])
    sA0 = sqrt(v.A0[end])
    bA0 = v.beta[end] / sA0

    fun(As) =
        As * v.R1 * (v.u[end] + sgamma * (ssAl - sqrt(sqrt(As)))) -
        (v.Pext + bA0 * (sqrt(As) - sA0)) + v.Pc

    dfun(As) = v.R1 * (v.u[end] + sgamma * (ssAl - 1.25 * sqrt(sqrt(As)))) - bA0 * 0.5 / sqrt(As)
    As = newtone(fun, dfun, As)
    us = (pressure(As, v.A0[end], v.beta[end], v.Pext) - v.Pout) / (As * v.R1)

    v.A[end] = As
    v.u[end] = us
end

function newtone(f::Function, df::Function, xn)
    for _=1:10
        xn -= f(xn) / df(xn)
    end
    xn
end

function update_ghost_cells!(n::Network)
    for v in values(n.vessels)
        v.U00A = v.A[1]
        v.U00Q = v.Q[1]
        v.UM1A = v.A[v.M]
        v.UM1Q = v.Q[v.M]
    end
end
