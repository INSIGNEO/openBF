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
    active_pane::Threads.Atomic{Int}   # 1=waveforms  2=convergence  3=log
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
                Threads.Atomic{Int}(1), Threads.Atomic{Int}(1), 0, 0, now, 0.0)
end

function record_step!(obs::TUIObserver, vessels, t::Float64, dt::Float64)
    obs.step_counter += 1
    obs.step_counter % obs.step_stride == 0 || return nothing
    @inbounds for i in eachindex(obs.monitored)
        spec = obs.monitored[i]
        push!(obs.waveforms[i].Q, Float32(vessels[spec.vessel_idx].Q[spec.node_idx]))
    end
    obs.scalars.t = t
    obs.scalars.dt = dt
    obs.scalars.gc_bytes = Base.gc_live_bytes()
    # cooperative yield so render task can run on single-threaded Julia.
    yield()
    nothing
end

function record_cycle!(obs::TUIObserver, cycle_idx::Int, vessels,
                       conv_P::Float64, conv_A::Float64 = 0.0, conv_Q::Float64 = 0.0)
    # guard against the floatmax() sentinel used before cycle 1
    conv_P > 0 && conv_P < floatmax() && push!(obs.convergence.P, Float32(conv_P))
    conv_A > 0 && push!(obs.convergence.A, Float32(conv_A))
    conv_Q > 0 && push!(obs.convergence.Q, Float32(conv_Q))
    obs.scalars.passed_cycles = cycle_idx
    # snapshot P waveforms from checkpoint data (col 4 = node3 = midpoint)
    # v.waveforms["P"] holds the just-completed cycle; v.waveforms_prev holds the one before
    @inbounds for i in eachindex(obs.monitored)
        v = vessels[obs.monitored[i].vessel_idx]
        haskey(v.waveforms, "P") || continue
        obs.waveforms[i].prev_P = cycle_idx > 0 ?
            Float32.(view(v.waveforms_prev["P"], :, 4)) .* _PA_TO_MMHG :
            Float32[]
        obs.waveforms[i].curr_P = Float32.(view(v.waveforms["P"], :, 4)) .* _PA_TO_MMHG
    end
    nothing
end
