const _WAVEFORM_SAMPLES = 400
const _PLOT_HEIGHT      = 6
const _PLOT_OVERHEAD    = 17  # y-axis label + borders + legend space

function _prewarm_plots()
    dummy = lineplot(zeros(Float64, 10); name="curr", height=_PLOT_HEIGHT, width=40,
                     canvas=BrailleCanvas, color=:cyan)
    lineplot!(dummy, zeros(Float64, 10); name="prev", color=:yellow)
    repr(MIME("text/plain"), dummy)
    nothing
end

function _waveform_panel(site::WaveformSite, panel_width::Int)
    p_live = latest(site.P, _WAVEFORM_SAMPLES)
    isempty(p_live) && return nothing

    plot_width = max(20, panel_width - _PLOT_OVERHEAD)
    prev       = site.prev_P

    # anchor y-axis to the full previous cycle so the live trace doesn't zoom in
    y_ref = isempty(prev) ? p_live : prev
    lo    = Float64(minimum(y_ref))
    hi    = Float64(maximum(y_ref))
    pad   = max((hi - lo) * 0.05, 1.0)  # ≥1 mmHg breathing room
    ylim  = (lo - pad, hi + pad)

    plt = lineplot(Float64.(p_live);
                   title  = "",
                   name   = "curr",
                   height = _PLOT_HEIGHT,
                   width  = plot_width,
                   ylim   = ylim,
                   canvas = BrailleCanvas,
                   color  = :cyan)

    if !isempty(prev)
        x_prev = collect(LinRange(1.0, Float64(length(p_live)), length(prev)))
        lineplot!(plt, x_prev, Float64.(prev); name = "prev", color = :yellow)
    end

    Panel(sprint(show, MIME("text/plain"), plt; context = :color => true);
          title   = site.label * "  P (mmHg)",
          width   = panel_width,
          padding = (0, 0, 0, 0),
          style   = "cyan dim")
end

function draw_waveforms(obs::TUIObserver, term_cols::Int = 80)
    isempty(obs.waveforms) && return nothing
    panel_width = max(40, term_cols - _SIDEBAR_WIDTH)
    idx = clamp(obs.active_site[], 1, length(obs.waveforms))
    _waveform_panel(obs.waveforms[idx], panel_width)
end
