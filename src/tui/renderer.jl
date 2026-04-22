const _RENDER_HZ = 30
const _RENDER_INTERVAL = 1.0 / _RENDER_HZ

function _render_frame(obs::TUIObserver)
    vals = latest(obs.waveforms[1].P, 80)
    isempty(vals) && return
    @printf "[TUI] t=%.4f  P̄=%6.1f  Pmin=%6.1f  Pmax=%6.1f\n" obs.scalars.t mean(vals) minimum(vals) maximum(vals)
end

function _render_loop(obs::TUIObserver)
    try
        while !obs.should_stop[]
            sleep(_RENDER_INTERVAL)
            obs.should_stop[] && break
            _render_frame(obs)
        end
    catch e
        e isa InterruptException || rethrow()
    end
    nothing
end

function start_render!(obs::TUIObserver)
    obs.should_stop[] = false
    obs.render_task = Threads.@spawn :interactive _render_loop(obs)
    nothing
end

function stop_render!(obs::TUIObserver)
    obs.should_stop[] = true
    t = obs.render_task
    t === nothing && return
    timedwait(() -> istaskdone(t), 2.0; pollint = 0.05)
    nothing
end
