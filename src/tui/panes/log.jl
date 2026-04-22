function draw_log(obs::TUIObserver, term_cols::Int, height::Int)
    panel_width = max(40, term_cols - _SIDEBAR_WIDTH)
    n_lines     = max(1, height - 2)  # subtract top/bottom borders
    entries     = latest(obs.log, n_lines)

    content = if isempty(entries)
        " (no log entries yet)"
    else
        join((@sprintf("%5.0fs  %s", t - obs.start_time, msg) for (t, msg) in entries), "\n")
    end

    Panel(content; title = "log", width = panel_width, height = height,
          padding = (0, 0, 0, 0), style = "dim")
end
