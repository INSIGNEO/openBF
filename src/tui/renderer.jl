const _RENDER_HZ       = 30
const _RENDER_INTERVAL = 1.0 / _RENDER_HZ
const _CURSOR_HIDE     = "\e[?25l"
const _CURSOR_SHOW     = "\e[?25h"
const _ALT_ENTER       = "\e[?1049h"
const _ALT_EXIT        = "\e[?1049l"
const _SCREEN_HOME     = "\e[H"
const _SCREEN_CLEAR    = "\e[2J\e[H"
const _BODY_HEIGHT     = _PLOT_HEIGHT + 4  # plot rows + border + title + legend

function _term_size()
    try displaysize(stdout) catch; (24, 80) end
end

function _update_rate!(obs::TUIObserver)
    now = time()
    elapsed = now - obs.last_check_time
    elapsed < 0.5 && return
    obs.steps_per_sec   = (obs.step_counter - obs.last_step_count) / elapsed
    obs.last_step_count = obs.step_counter
    obs.last_check_time = now
end

const _PANE_LABELS = ("waveforms", "convergence", "log")

function _render_frame(obs::TUIObserver, buf::IOBuffer)
    _update_rate!(obs)
    obs.frame_count += 1
    rows, cols = _term_size()

    header  = _make_header(obs, cols)
    sidebar = _make_sidebar(obs, _BODY_HEIGHT)
    pane    = clamp(obs.active_pane[], 1, length(_PANE_LABELS))

    body = if pane == 2
        draw_convergence(obs, cols)
    elseif pane == 3
        draw_log(obs, cols, _BODY_HEIGHT)
    else
        draw_waveforms(obs, cols)
    end

    truncate(buf, 0); seek(buf, 0)
    ctx = IOContext(buf, :color => true, :displaysize => (rows, cols))
    print(ctx, _SCREEN_HOME)
    if body !== nothing
        print(ctx, header / (body * sidebar))
    else
        label = pane == 2 ?
            " Waiting for convergence data… (need ≥1 completed cycle)" :
            " Waiting for simulation data… (stride=$(obs.step_stride) steps)"
        waiting = Panel(label; width = cols, padding = (0,0,0,0), style = "dim")
        print(ctx, header / waiting)
    end
    write(stdout, take!(buf))
    flush(stdout)
    nothing
end

function _render_loop(obs::TUIObserver)
    buf     = IOBuffer()
    err_msg = ""
    print(stdout, _SCREEN_CLEAR); flush(stdout)
    try
        while !obs.should_stop[]
            sleep(_RENDER_INTERVAL)
            obs.should_stop[] && break
            try
                _render_frame(obs, buf)
                err_msg = ""  # clear any previous error once a frame succeeds
            catch e
                msg = sprint(showerror, e)
                if msg != err_msg
                    err_msg = msg
                    truncate(buf, 0); seek(buf, 0)
                    write(buf, _SCREEN_HOME)
                    write(buf, _SCREEN_CLEAR)
                    write(buf, "TUI render error (solver still running):\n")
                    write(buf, msg)
                    write(buf, "\n")
                    write(stdout, take!(buf)); flush(stdout)
                end
            end
        end
    catch e
        e isa InterruptException || rethrow()
    end
    nothing
end

function start_render!(obs::TUIObserver)
    _prewarm_plots()
    install_log_sink!(obs)
    print(stdout, _ALT_ENTER, _CURSOR_HIDE); flush(stdout)
    obs.should_stop[] = false
    obs.render_task = if Threads.nthreads(:interactive) > 0
        Threads.@spawn :interactive _render_loop(obs)
    else
        Threads.@spawn :default _render_loop(obs)
    end
    nothing
end

function stop_render!(obs::TUIObserver)
    obs.should_stop[] = true
    t = obs.render_task
    t !== nothing && timedwait(() -> istaskdone(t), 2.0; pollint = 0.05)
    print(stdout, _CURSOR_SHOW, _ALT_EXIT); flush(stdout)
    uninstall_log_sink!(obs)
    nothing
end
