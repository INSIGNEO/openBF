const _WAVEFORM_SAMPLES = 400
const _PLOT_HEIGHT = 6

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
          width   = plot_width + 12,
          padding = (0, 0, 0, 0),
          style   = "cyan dim")
end

function draw_waveforms(obs::TUIObserver, io::IO, term_cols::Int = 80)
    plot_width = max(20, term_cols - 12)
    panels = []
    for site in obs.waveforms
        p = _waveform_panel(site, plot_width)
        p === nothing && continue
        push!(panels, p)
    end
    isempty(panels) && return false
    stacked = reduce(/, panels)
    print(io, stacked)
    true
end
