const _SIDEBAR_WIDTH = 26
const _PA_TO_MMHG    = 1.0 / 133.322

function _make_sidebar(obs::TUIObserver, height::Int)
    s      = obs.scalars
    idx    = clamp(obs.active_site[], 1, length(obs.waveforms))
    site   = obs.waveforms[idx]
    p_vals = latest(site.P, 400)
    q_vals = latest(site.Q, 400)

    map_p  = isempty(p_vals) ? 0.0 : Float64(mean(p_vals))
    peak_p = isempty(p_vals) ? 0.0 : Float64(maximum(p_vals))
    pp     = isempty(p_vals) ? 0.0 : Float64(maximum(p_vals) - minimum(p_vals))
    q_m    = isempty(q_vals) ? 0.0 : Float64(mean(q_vals))

    lines = (
        @sprintf("MAP  %7.1f mmHg", map_p  * _PA_TO_MMHG),
        @sprintf("PP   %7.1f mmHg", pp     * _PA_TO_MMHG),
        @sprintf("Peak %7.1f mmHg", peak_p * _PA_TO_MMHG),
        @sprintf("Qmn  %8.5f m3/s", q_m),
        @sprintf("dt   %8.2e s",    s.dt),
        @sprintf("GC   %7.1f MB",   s.gc_bytes / 1e6),
        @sprintf("cyc  %7d",        s.passed_cycles),
    )

    Panel(join(lines, "\n"); title = "vitals", width = _SIDEBAR_WIDTH,
          height = height, padding = (0,0,0,0), style = "green dim")
end
