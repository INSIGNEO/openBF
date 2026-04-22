const _CONV_SAMPLES = 512

function draw_convergence(obs::TUIObserver, term_cols::Int = 80)
    p_vals = latest(obs.convergence.P, _CONV_SAMPLES)
    isempty(p_vals) && return nothing

    panel_width = max(40, term_cols - _SIDEBAR_WIDTH)
    plot_width  = max(20, panel_width - _PLOT_OVERHEAD)

    # log scale requires strictly positive values
    y_p = clamp.(Float64.(p_vals), 1e-10, Inf)
    plt = lineplot(y_p;
                   name   = "P",
                   height = _PLOT_HEIGHT,
                   width  = plot_width,
                   yscale = :log10,
                   canvas = BrailleCanvas,
                   color  = :cyan)

    q_vals = latest(obs.convergence.Q, _CONV_SAMPLES)
    if !isempty(q_vals) && any(>(0f0), q_vals)
        lineplot!(plt, clamp.(Float64.(q_vals), 1e-10, Inf); name = "Q", color = :magenta)
    end

    a_vals = latest(obs.convergence.A, _CONV_SAMPLES)
    if !isempty(a_vals) && any(>(0f0), a_vals)
        lineplot!(plt, clamp.(Float64.(a_vals), 1e-10, Inf); name = "A", color = :yellow)
    end

    Panel(sprint(show, MIME("text/plain"), plt; context = :color => true);
          title   = "convergence  RMSE/cycle",
          width   = panel_width,
          padding = (0, 0, 0, 0),
          style   = "cyan dim")
end
