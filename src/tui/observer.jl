record_step!(::Nothing, args...) = nothing
record_cycle!(::Nothing, args...) = nothing
record_log!(::Nothing, args...) = nothing
start_render!(::Nothing) = nothing
stop_render!(::Nothing) = nothing

# ---------------------------------------------------------------------------

struct MonitorSpec
    vessel_idx::Int
    node_idx::Int
    label::String
end

mutable struct TUIObserver
    waveforms::Vector{WaveformSite}
    convergence::ConvergenceHistory
    scalars::Snapshot
    log::LogRing
    monitored::Vector{MonitorSpec}
    step_stride::Int
    step_counter::Int
    should_stop::Threads.Atomic{Bool}
    render_task::Union{Task, Nothing}
    # sim metadata (set at construction, read-only after)
    sim_name::String
    t_final::Float64
    n_cycles::Int
    start_time::Float64
    # render-thread state (only written by render loop)
    active_site::Threads.Atomic{Int}
    frame_count::Int
    last_step_count::Int
    last_check_time::Float64
    steps_per_sec::Float64
end

const _PRIORITY_PATTERNS = ("aorta", "carotid", "brachial", "femoral")

function _select_monitored(vessels_vec, is_outlet)
    matched = MonitorSpec[]
    for pat in _PRIORITY_PATTERNS
        for (i, v) in enumerate(vessels_vec)
            occursin(pat, lowercase(String(v.label))) || continue
            any(s -> s.vessel_idx == i, matched) && continue
            push!(matched, MonitorSpec(i, v.node3, String(v.label)))
            length(matched) >= 4 && return matched
            break
        end
    end
    isempty(matched) || return matched
    # fallback: first 4 outlet vessels
    for (i, v) in enumerate(vessels_vec)
        is_outlet[i] || continue
        push!(matched, MonitorSpec(i, v.node3, String(v.label)))
        length(matched) >= 4 && return matched
    end
    isempty(matched) || return matched
    # absolute fallback: first 4 vessels
    for i in 1:min(4, length(vessels_vec))
        v = vessels_vec[i]
        push!(matched, MonitorSpec(i, v.node3, String(v.label)))
    end
    matched
end

function TUIObserver(vessels_vec, is_outlet;
                     step_stride::Int = 200,
                     sim_name::String = "",
                     t_final::Float64 = 0.0,
                     n_cycles::Int = 0)
    monitored = _select_monitored(vessels_vec, is_outlet)
    waveforms = [WaveformSite(spec.label) for spec in monitored]
    now = time()
    TUIObserver(waveforms, ConvergenceHistory(), Snapshot(), LogRing(128), monitored,
                step_stride, 0, Threads.Atomic{Bool}(false), nothing,
                sim_name, t_final, n_cycles, now,
                Threads.Atomic{Int}(1), 0, 0, now, 0.0)
end

function record_step!(obs::TUIObserver, vessels, t::Float64, dt::Float64)
    obs.step_counter += 1
    obs.step_counter % obs.step_stride == 0 || return nothing
    @inbounds for i in eachindex(obs.monitored)
        spec = obs.monitored[i]
        v = vessels[spec.vessel_idx]
        push!(obs.waveforms[i].P, Float32(v.P[spec.node_idx]))
        push!(obs.waveforms[i].Q, Float32(v.Q[spec.node_idx]))
    end
    obs.scalars.t = t
    obs.scalars.dt = dt
    obs.scalars.gc_bytes = Base.gc_live_bytes()
    nothing
end

function record_cycle!(obs::TUIObserver, cycle_idx::Int, conv_error::Float64)
    push!(obs.convergence.P, Float32(conv_error))
    obs.scalars.passed_cycles = cycle_idx
    nothing
end
