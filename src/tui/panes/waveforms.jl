const _WAVEFORM_SAMPLES = 400
const _PLOT_HEIGHT      = 6

function _prewarm_plots()
    dummy = lineplot(zeros(Float64, 10); name="P", height=_PLOT_HEIGHT, width=40,
                     canvas=BrailleCanvas, color=:cyan)
    lineplot!(dummy, zeros(Float64, 10); name="Q", color=:magenta)
    repr(MIME("text/plain"), dummy)
    nothing
end

function _waveform_panel(site::WaveformSite, plot_width::Int)
    p_vals = latest(site.P, _WAVEFORM_SAMPLES)
    isempty(p_vals) && return nothing
    q_vals = latest(site.Q, _WAVEFORM_SAMPLES)

    plt = lineplot(Float64.(p_vals);
                   title  = "",
                   name   = "P",
                   height = _PLOT_HEIGHT,
                   width  = plot_width,
                   canvas = BrailleCanvas,
                   color  = :cyan)
    if length(q_vals) == length(p_vals)
        lineplot!(plt, Float64.(q_vals); name = "Q", color = :magenta)
    end

    Panel(repr(MIME("text/plain"), plt);
          title   = site.label,
          width   = plot_width,
          padding = (0, 0, 0, 0),
          style   = "cyan dim")
end

# Returns the Panel for the active site only (or nothing if no data yet).
function draw_waveforms(obs::TUIObserver, term_cols::Int = 80)
    isempty(obs.waveforms) && return nothing
    plot_width = max(20, term_cols - _SIDEBAR_WIDTH - 2)
    idx  = clamp(obs.active_site[], 1, length(obs.waveforms))
    _waveform_panel(obs.waveforms[idx], plot_width)
end
