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

Parameterised in `(u, α)` where `α = A^{1/4}`, matching the three legacy
solvers exactly. Unknowns are laid out as `[u₁…uₖ, α₁…αₖ]`.

# Fields
- `id::Int`: unique identifier (the shared node ID) for diagnostics.
- `vessels::Vector{Int}`: indices into the network's vessel list.
- `sides::Vector{Symbol}`: `:outlet` or `:inlet`, per attached vessel.
- `signs::Vector{Int}`: +1 for `:outlet`, −1 for `:inlet`. Precomputed.
- `use_total_pressure::Bool`: when true, the pressure-continuity equation
  uses total pressure (static + dynamic); set to `true` for k=2 (conjunction).

Work arrays (preallocated to avoid per-step allocation):

- `x::Vector{Float64}`: unknowns, length 2k, ordered [u₁…uₖ, α₁…αₖ].
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

# Returns (u, α) where α = A^{1/4} from the vessel's boundary face.
face_u_alpha(v::Vessel, side::Symbol) =
    side == :outlet ? (v.u[end], v.A[end]^0.25) : (v.u[1], v.A[1]^0.25)

# Wave-speed coefficient: c = k·α  (matches legacy sqrt(1.5·gamma))
face_k(v::Vessel, side::Symbol) =
    side == :outlet ? sqrt(1.5 * v.gamma[end]) : sqrt(1.5 * v.gamma[1])

face_A0_beta(v::Vessel, side::Symbol) =
    side == :outlet ? (v.A0[end], v.beta[end]) : (v.A0[1], v.beta[1])

# ── Newton loop ───────────────────────────────────────────────────────────────

function pack_initial_guess!(jc::Junction, vessels::Vector{Vessel})
    k = length(jc.vessels)
    for (i, vid) in enumerate(jc.vessels)
        v = vessels[vid]
        u, alpha = face_u_alpha(v, jc.sides[i])
        jc.x[i]   = u
        jc.x[k+i] = alpha
    end
    return jc
end

function residual!(jc::Junction, vessels::Vector{Vessel}, rho::Float64,
                   Wstar::Vector{Float64})
    k = length(jc.vessels)

    # Block 1 — characteristic equations (linear in α, matching legacy)
    for i in 1:k
        v    = vessels[jc.vessels[i]]
        ki   = face_k(v, jc.sides[i])
        jc.F[i] = jc.x[i] + jc.signs[i] * 4ki * jc.x[k+i] - Wstar[i]
    end

    # Block 2 — mass conservation: Σ sᵢ · uᵢ · αᵢ⁴ = 0
    mass = 0.0
    for i in 1:k
        mass += jc.signs[i] * jc.x[i] * jc.x[k+i]^4
    end
    jc.F[k+1] = mass

    # Block 3 — pressure continuity: Pⱼ - P₁ = 0  (j = 2…k)
    v1        = vessels[jc.vessels[1]]
    a1        = jc.x[k+1]
    A01, β1   = face_A0_beta(v1, jc.sides[1])
    P1        = β1 * (a1^2 / sqrt(A01) - 1.0)
    jc.use_total_pressure && (P1 += 0.5rho * jc.x[1]^2)

    for j in 1:k-1
        vj       = vessels[jc.vessels[j+1]]
        aj       = jc.x[k+j+1]
        A0j, βj  = face_A0_beta(vj, jc.sides[j+1])
        Pj       = βj * (aj^2 / sqrt(A0j) - 1.0)
        jc.use_total_pressure && (Pj += 0.5rho * jc.x[j+1]^2)
        jc.F[k+1+j] = Pj - P1
    end

    return jc.F
end

function jacobian!(jc::Junction, vessels::Vector{Vessel}, rho::Float64)
    k = length(jc.vessels)
    fill!(jc.J, 0.0)

    # Characteristic rows — diagonal in (uᵢ, αᵢ)
    for i in 1:k
        v  = vessels[jc.vessels[i]]
        ki = face_k(v, jc.sides[i])
        jc.J[i, i]   = 1.0
        jc.J[i, k+i] = jc.signs[i] * 4ki
    end

    # Mass conservation row
    for i in 1:k
        u_i = jc.x[i]
        a_i = jc.x[k+i]
        jc.J[k+1, i]   = jc.signs[i] * a_i^4
        jc.J[k+1, k+i] = jc.signs[i] * 4u_i * a_i^3
    end

    # Pressure continuity rows
    v1        = vessels[jc.vessels[1]]
    a1        = jc.x[k+1]
    u1        = jc.x[1]
    A01, β1   = face_A0_beta(v1, jc.sides[1])
    dP1da     = 2β1 * a1 / sqrt(A01)

    for j in 1:k-1
        vj       = vessels[jc.vessels[j+1]]
        aj       = jc.x[k+j+1]
        uj       = jc.x[j+1]
        A0j, βj  = face_A0_beta(vj, jc.sides[j+1])
        dPjda    = 2βj * aj / sqrt(A0j)
        row      = k + 1 + j

        jc.J[row, k+j+1] =  dPjda
        jc.J[row, k+1]   = -dP1da

        if jc.use_total_pressure
            jc.J[row, j+1] =  rho * uj
            jc.J[row, 1]   = -rho * u1
        end
    end

    return jc.J
end

function unpack_solution!(jc::Junction, vessels::Vector{Vessel})
    k = length(jc.vessels)
    for (i, vid) in enumerate(jc.vessels)
        v     = vessels[vid]
        u     = jc.x[i]
        alpha = jc.x[k+i]
        A     = alpha^4
        Q     = u * A
        if jc.sides[i] == :outlet
            v.u[end] = u;  v.A[end] = A;  v.Q[end] = Q
        else
            v.u[1]   = u;  v.A[1]   = A;  v.Q[1]   = Q
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
    k = length(jc.vessels)
    pack_initial_guess!(jc, vessels)

    # Precompute Riemann invariants from the current (pre-Newton) face state.
    Wstar = ntuple(k) do i
        v  = vessels[jc.vessels[i]]
        ki = face_k(v, jc.sides[i])
        jc.x[i] + jc.signs[i] * 4ki * jc.x[k+i]
    end
    Wstar_vec = collect(Float64, Wstar)

    converged = false
    for _ in 1:JUNCTION_NEWTON_MAXIT
        residual!(jc, vessels, blood.rho, Wstar_vec)
        norm(jc.F) < JUNCTION_NEWTON_TOL && (converged = true; break)
        jacobian!(jc, vessels, blood.rho)
        jc.dx .= jc.J \ jc.F
        jc.x  .-= jc.dx
    end
    converged || @warn "junction $(jc.id) did not converge" norm(jc.F)
    unpack_solution!(jc, vessels)
    return jc
end
