# AGENTS.md

Guidance for AI coding agents working on openBF. Human contributors are welcome to read this too, but the intended audience is agents.

---

## What openBF is

openBF is a 1D hemodynamics solver. It simulates blood flow through networks of arbitrary topology using a MUSCL finite-volume scheme on a directed graph of 1D vessels. Junctions (bifurcations, conjunctions, anastomoses) are resolved with Newton iteration on small (4×4, 6×6) systems. Outlets use either a three-element Windkessel model or a characteristic-based compatibility condition.

The math is correct. The software engineering has room to move. Assume any change you propose is a performance or clarity change, not a physics change.

---

## Repository layout

```
src/
  openBF.jl              module entry
  vessel.jl              Vessel struct + constructors
  network.jl             Network + Heart
  solver.jl              MUSCL time-stepping (hot path)
  boundary_conditions.jl inlet, outlet, WK3
  bifurcations.jl        3-vessel junction Newton solve
  conjunctions.jl        2-vessel junction Newton solve
  anastomosis.jl         inverse 3-vessel junction Newton solve
  simulation.jl          outer loop, convergence, I/O orchestration
  output.jl              waveform saving
models/                  YAML network configs + inlet waveforms
test/                    correctness tests
bench/                   benchmarking harness (see §Benchmarking)
```

Read in that order before making your first change. Do not skip this; the mental model matters more than the specific bug you're chasing.

---

## Orientation checklist

Before your first edit, confirm you understand:

- The graph topology does not change mid-simulation. Treat it as compile-time data.
- `Vessel` holds both immutable geometry (`M`, `beta`, `A0`, `gamma`, …) and mutable state (`A`, `Q`, `u`, `P`).
- `solve!` in `solver.jl` is called once per time step and is where hot-path work happens.
- Junctions are solved by iterating over edges, not nodes. A bifurcation solve writes across three vessels simultaneously.
- `run_simulation` in `simulation.jl` does real work: it `cd`s into a results directory, reads YAML, constructs the network, and loops until convergence or a cycle cap.

If any of these statements surprises you, re-read the relevant file.

---

## Core principles

These are non-negotiable. If a proposed change violates one, find another way.

### KISS

Smallest diff that moves the needle. No speculative abstractions. No rewrites when a patch suffices. No trait systems, no generated code, no metaprogramming unless the alternative is materially worse.

### DRY, but earned

When you see duplicated code, factor it *once*, only when factoring reduces line count or fixes a real maintenance burden. Do not factor speculatively "for future flexibility". Duplicated code with different invariants is not duplication.

### Atomic commits

One idea per commit. One measurable effect. If a diff is growing past ~80 lines or touching more than 3 files, split it. If a step is titled "X and also Y", it's two commits.

### Benchmark-gated

Every performance change must be validated against the harness in §Benchmarking *before* moving to the next change. A change that doesn't beat the noise floor gets reverted or merged with explicit justification (e.g., "enables next step"). No unmeasured optimizations.

### Correctness-gated

Every change must pass the regression test in §Correctness. Default tolerance is bit-exact. Steps that relax tolerance must state so explicitly and update the reference.

### When in doubt, don't

If a change's benefit is unclear, the diff is getting ugly, or you're three layers deep in a refactor that wasn't in your plan, stop. Report what you found. Ask. Do not bulldoze.

---

## Idiomatic Julia conventions

- `snake_case` for functions and variables; `CamelCase` for types and modules.
- Prefer `struct` over `mutable struct`. Use `mutable struct` only for state that genuinely evolves in place.
- Concrete types on struct fields in any type that lives in the hot path. No `Dict()` without type parameters.
- `@view` for array slices in hot code. `A[:, j]` allocates; `@view A[:, j]` does not.
- `eachindex(x)` over `1:length(x)`. `eachindex(x, y)` when iterating two aligned arrays.
- `@inbounds` only after you've verified bounds are safe. `@simd` only on loops with no loop-carried dependencies and simple bodies.
- Parametric functions over `::Function` arguments. `f(g::G) where {G}` specializes; `f(g::Function)` boxes.
- `SVector`/`SMatrix` for small fixed-size linear algebra. Avoid the `zeros(6,6)` → `SMatrix(...)` anti-pattern — build `MMatrix` then convert, or construct directly from a tuple.
- Prefer `ntuple(f, Val(N))` over `ntuple(f, N)` when `N` is a compile-time constant.
- Return `nothing` explicitly from mutating functions (`foo!`). Use `return` at the bottom, not a bare expression.
- `const` for module-level constants. No magic numbers in hot loops — bind them to names.
- Docstrings on exported functions. No comments explaining *what* the code does; comments should explain *why* when non-obvious.

When in doubt, match the style of the surrounding code. If the surrounding code is bad, you may fix it — but in a separate commit.

---

## Benchmarking

### Harness

Benchmarks live in `bench/harness.jl`. The three canonical models are:

- `cca` — single vessel. Proves MUSCL kernel changes.
- `ibif` — one bifurcation. Proves junction changes.
- `adan56` — 56-vessel arterial network. Proves network/dispatch/threading changes.

A change that only helps `cca` is a kernel win. A change that only helps `adan56` is a network win. Ideally both.

### Workflow

For every performance change:

1. Before the change: `run_suite("before_<name>")`.
2. Apply change.
3. After: `run_suite("after_<name>")`.
4. Compare: `judge_against("after_<name>", "before_<name>")`.
5. Record the result in the commit message.

### Noise floor

On first setup, and after any environment change, run `run_suite` twice back-to-back without code changes. The |Δ%| between runs is your noise floor — typically 1–3% on a quiet machine. A change must beat this by **at least 2×** to count as real. Smaller improvements may still be merged if they're structurally motivated, but must be flagged as such.

### Regression policy

A change must:
- Improve at least one benchmark by ≥2× the noise floor, OR structurally enable a future change.
- Not regress any benchmark by more than 1× the noise floor.

Otherwise: revert or reformulate.

### Threading

Any change touching parallelism or shared state must be benchmarked with both `JULIA_NUM_THREADS=1` and `JULIA_NUM_THREADS=$(nproc)`. Single-threaded performance must not regress.

### Profiling

Every 3–4 changes, confirm the hot path has shifted:

```julia
using Profile, ProfileView
bench_model("adan56")  # warm up
Profile.clear(); @profile bench_model("adan56"); ProfileView.view()
```

Record the top-3 functions by self-time in the commit message of the next change. If the hot path hasn't shifted after several changes, you're working on the wrong thing.

---

## Correctness

### Reference capture

`test/reference.jl` captures the final `A`, `Q`, `u`, `P` arrays of every vessel across all three canonical models after a short, deterministic run (3 cycles, convergence tolerance disabled). References are serialized to `test/ref_<model>.jls`.

Regenerate references only when a change explicitly alters behavior (e.g., fixing a bug in a formula). Document the reason in the commit message.

### Regression test

```julia
check_regression("cca"; rtol=0.0, atol=0.0)     # bit-exact default
check_regression("ibif"; rtol=1e-10)            # after FP-reordering changes
```

### Tolerance policy

- **Refactors and structural changes: bit-exact** (`rtol=0, atol=0`). If it's not bit-exact, you have a bug or a hidden semantic change.
- **`@fastmath` or reassociation changes:** `rtol=1e-8`, scoped narrowly, verified per-field.
- **Intentional behavior changes:** no tolerance check; regenerate reference and document.
- **End of a phase:** may loosen default to `rtol=1e-10` to account for accumulated FP reordering across many small changes.

### Known correctness issues in the current code

These are pre-existing bugs. Do not introduce more.

- `check()` in `network.jl` has a `||`/`&&` precedence trap that silences the "too many edges at vertex" error in some cases.
- `get_conv_error` computes `|mean(error)|`, not RMSE, despite its name and the "RMSE (mmHg)" label.

If you find more, flag them. Do not silently fix — file an issue first, fix in a dedicated commit with clear tolerance discussion.

---

## What NOT to do

- Do not introduce new dependencies to `Project.toml` without strong justification. `JSON` for the benchmark harness is fine; adding a new linear algebra library is not.
- Do not rewrite working code for aesthetics.
- Do not refactor files you aren't actively changing.
- Do not add comments explaining what code does — explain *why*, only when non-obvious.
- Do not rename things "while you're here". Separate commit.
- Do not change physical formulas (Olufsen wall-thickness coefficients, MUSCL flux, WK3 equations, Newton system definitions). If one looks wrong, file an issue and stop.
- Do not touch the YAML schema or config parsing. Users have configs in the wild.
- Do not GPU-port or add threading without an explicit plan that survived benchmark-gating.
- Do not merge a change whose benchmark data you haven't personally captured.

---

## Stop conditions

Halt and report — do not continue — if:

- A correctness regression appears that you cannot resolve in ~30 minutes of focused work.
- A change regresses performance on any benchmark beyond the noise floor and one alternative formulation also failed.
- A diff is growing past ~150 lines despite your best efforts to keep it atomic.
- You find yourself mid-refactor doing something not in your plan.
- You hit an ambiguity in the physics or the numerical scheme that isn't resolved by reading the code.

"Report" means: summarize the state, the data, the proposed alternatives, and ask.

---

## Commit message format

```
<phase/step or area>: <one-line summary>

<what changed and why, 2–5 lines>

Benchmarks (vs. previous commit, min time, n=5):
  cca:    -X.X%  (Yk→Zk allocs)
  ibif:   -X.X%
  adan56: -X.X%

Correctness: bit-exact | rtol=1e-N | reference regenerated (reason)
Profile top-3 (if re-profiled): fn1 (X%), fn2 (Y%), fn3 (Z%)
```

Empty benchmark sections are a smell. Either you didn't run them, or the change wasn't benchmark-relevant — say which.

---

## Final deliverable (for larger efforts)

When completing a multi-step plan, provide:

1. A table: `step | cca_Δ% | ibif_Δ% | adan56_Δ% | alloc_Δ% | notes`.
2. Cumulative speedup vs. the starting commit, all three benchmarks.
3. Skipped steps with one-line rationale each.
4. Any correctness issues found (in the original code or introduced and fixed mid-stream).
5. Updated reference tolerances and rationale.

---

## When this file is wrong

This file encodes current beliefs. If you discover something in it is wrong, outdated, or actively misleading, update it — in a dedicated commit, separate from your other work. Keep the change minimal and the diff readable.
