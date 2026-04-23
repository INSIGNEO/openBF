#=
Copyright 2015-2026 INSIGNEO Institute for in silico Medicine

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

# NOTE: Junction holds preallocated work arrays (x, F, J, dx). It is NOT
# safe to solve the same Junction from multiple threads concurrently.
# If the junction loop is ever parallelised, each thread must own a
# private Junction copy (or these arrays must move to a thread-local
# workspace).

"""
    Junction

A generic n-furcation: k ≥ 2 vessels meeting at a single point with
arbitrary parent (`:outlet`) / daughter (`:inlet`) split.

# Fields
- `id::Int`: unique identifier (the shared node ID) for diagnostics.
- `vessels::Vector{Int}`: indices into the network's vessel list.
- `sides::Vector{Symbol}`: `:outlet` or `:inlet`, per attached vessel.
- `signs::Vector{Int}`: +1 for `:outlet`, −1 for `:inlet`. Precomputed.
- `use_total_pressure::Bool`: when true, the pressure-continuity equation
  uses total pressure (static + dynamic); default `false` matches the
  bifurcation and anastomosis solvers. Set to `true` to match the
  conjunction solver.

Work arrays (preallocated to avoid per-step allocation):

- `x::Vector{Float64}`: unknowns, length 2k, ordered [A₁,Q₁,A₂,Q₂,…].
- `F::Vector{Float64}`: residual, length 2k.
- `J::Matrix{Float64}`: Jacobian, 2k × 2k.
- `dx::Vector{Float64}`: Newton update, length 2k.
"""
struct Junction
    id::Int
    vessels::Vector{Int}
    sides::Vector{Symbol}
    signs::Vector{Int}
    use_total_pressure::Bool
    x::Vector{Float64}
    F::Vector{Float64}
    J::Matrix{Float64}
    dx::Vector{Float64}
end

function Junction(id::Int, vessels::Vector{Int}, sides::Vector{Symbol};
                  use_total_pressure::Bool=false)
    length(vessels) == length(sides) || throw(ArgumentError(
        "vessels and sides must have equal length"))
    length(vessels) >= 2 || throw(ArgumentError(
        "a junction needs at least 2 attached vessels"))
    all(s -> s in (:inlet, :outlet), sides) || throw(ArgumentError(
        "sides must be :inlet or :outlet"))
    @assert any(s == :outlet for s in sides) "junction $(id): no outlet vessel"
    @assert any(s == :inlet  for s in sides) "junction $(id): no inlet vessel"
    k = length(vessels)
    signs = [s == :outlet ? +1 : -1 for s in sides]
    Junction(id, vessels, sides, signs, use_total_pressure,
             zeros(2k), zeros(2k), zeros(2k, 2k), zeros(2k))
end

# ── private helpers ───────────────────────────────────────────────────────────

face_state(v::Vessel, side::Symbol) =
    side == :outlet ? (v.A[end], v.Q[end]) : (v.A[1], v.Q[1])

face_gamma(v::Vessel, side::Symbol) =
    side == :outlet ? v.gamma[end] : v.gamma[1]

face_A0_beta(v::Vessel, side::Symbol) =
    side == :outlet ? (v.A0[end], v.beta[end]) : (v.A0[1], v.beta[1])

# Riemann invariant extrapolated from the interior to the face.
# Fixed before Newton starts; computed from the vessel's current boundary cell.
function extract_Wstar(v::Vessel, side::Symbol)
    A, Q  = face_state(v, side)
    gamma = face_gamma(v, side)
    u     = Q / A
    c     = wave_speed(A, gamma)
    side == :outlet ? u + 4c : u - 4c
end

# Pressure without Pext (Pext cancels in the continuity equation).
_static_P(A::Float64, A0::Float64, beta::Float64) = beta * (sqrt(A / A0) - 1.0)

function _junction_P(A::Float64, Q::Float64, v::Vessel, side::Symbol,
                     rho::Float64, use_total::Bool)
    A0, beta = face_A0_beta(v, side)
    p = _static_P(A, A0, beta)
    use_total ? p + 0.5 * rho * (Q / A)^2 : p
end

# ── Newton loop ───────────────────────────────────────────────────────────────

function pack_initial_guess!(jc::Junction, vessels::Vector{Vessel})
    for (i, vid) in enumerate(jc.vessels)
        v = vessels[vid]
        A, Q = face_state(v, jc.sides[i])
        jc.x[2i-1] = A
        jc.x[2i]   = Q
    end
    return jc
end

function residual!(jc::Junction, vessels::Vector{Vessel}, rho::Float64)
    k = length(jc.vessels)
    # Block 1 — characteristic equations, one per vessel
    for i in 1:k
        v     = vessels[jc.vessels[i]]
        A     = jc.x[2i-1]
        Q     = jc.x[2i]
        gamma = face_gamma(v, jc.sides[i])
        c     = wave_speed(A, gamma)
        Wstar = extract_Wstar(v, jc.sides[i])
        u     = Q / A
        jc.F[i] = jc.sides[i] == :outlet ? u + 4c - Wstar : u - 4c - Wstar
    end
    # Block 2 — mass conservation
    mass = 0.0
    for i in 1:k
        mass += jc.signs[i] * jc.x[2i]
    end
    jc.F[k+1] = mass
    # Block 3 — pressure continuity
    v1 = vessels[jc.vessels[1]]
    P1 = _junction_P(jc.x[1], jc.x[2], v1, jc.sides[1], rho, jc.use_total_pressure)
    for j in 1:k-1
        vj = vessels[jc.vessels[j+1]]
        Pj = _junction_P(jc.x[2(j+1)-1], jc.x[2(j+1)], vj, jc.sides[j+1], rho,
                         jc.use_total_pressure)
        jc.F[k+1+j] = Pj - P1
    end
    return jc.F
end

function jacobian!(jc::Junction, vessels::Vector{Vessel}, rho::Float64)
    k = length(jc.vessels)
    fill!(jc.J, 0.0)
    # Characteristic rows — couple only (Aᵢ, Qᵢ)
    for i in 1:k
        v     = vessels[jc.vessels[i]]
        A     = jc.x[2i-1]
        Q     = jc.x[2i]
        gamma = face_gamma(v, jc.sides[i])
        c     = wave_speed(A, gamma)
        sigma = Float64(jc.signs[i])
        jc.J[i, 2i-1] = -Q / A^2 + sigma * c / A   # ∂(u ± 4c)/∂A
        jc.J[i, 2i]   = 1.0 / A                      # ∂(u ± 4c)/∂Q
    end
    # Mass conservation row — couples the Qᵢ only
    for i in 1:k
        jc.J[k+1, 2i] = Float64(jc.signs[i])
    end
    # Pressure continuity rows — couple (Aⱼ₊₁, Qⱼ₊₁) and (A₁, Q₁)
    v1        = vessels[jc.vessels[1]]
    A1, Q1    = jc.x[1], jc.x[2]
    A01, beta1 = face_A0_beta(v1, jc.sides[1])
    dP1dA     = beta1 / (2sqrt(A1 * A01))
    for j in 1:k-1
        vj         = vessels[jc.vessels[j+1]]
        Aj         = jc.x[2(j+1)-1]
        Qj         = jc.x[2(j+1)]
        A0j, betaj = face_A0_beta(vj, jc.sides[j+1])
        dPjdA      = betaj / (2sqrt(Aj * A0j))
        row        = k + 1 + j
        if jc.use_total_pressure
            jc.J[row, 2(j+1)-1] =  dPjdA - rho * Qj^2 / Aj^3
            jc.J[row, 2(j+1)]   =  rho * Qj / Aj^2
            jc.J[row, 1]        = -(dP1dA - rho * Q1^2 / A1^3)
            jc.J[row, 2]        = -rho * Q1 / A1^2
        else
            jc.J[row, 2(j+1)-1] =  dPjdA
            jc.J[row, 1]        = -dP1dA
        end
    end
    return jc.J
end

function unpack_solution!(jc::Junction, vessels::Vector{Vessel})
    for (i, vid) in enumerate(jc.vessels)
        v = vessels[vid]
        A = jc.x[2i-1]
        Q = jc.x[2i]
        if jc.sides[i] == :outlet
            v.A[end] = A
            v.Q[end] = Q
            v.u[end] = Q / A
        else
            v.A[1] = A
            v.Q[1] = Q
            v.u[1] = Q / A
        end
    end
end

const JUNCTION_NEWTON_TOL   = 1e-5   # matches NRbif / NRconj / NRan tolerance
const JUNCTION_NEWTON_MAXIT = 30     # matches legacy maxiter

"""
    solve_junction!(jc, vessels, blood)

Solve the n-furcation junction `jc` in-place using Newton's method and
write the converged face state back into each attached vessel.
"""
function solve_junction!(jc::Junction, vessels::Vector{Vessel}, blood::Blood)
    pack_initial_guess!(jc, vessels)
    converged = false
    for _ in 1:JUNCTION_NEWTON_MAXIT
        residual!(jc, vessels, blood.rho)
        norm(jc.F) < JUNCTION_NEWTON_TOL && (converged = true; break)
        jacobian!(jc, vessels, blood.rho)
        jc.dx .= jc.J \ jc.F
        jc.x  .-= jc.dx
    end
    converged || @warn "junction $(jc.id) did not converge" norm(jc.F)
    unpack_solution!(jc, vessels)
    return jc
end
