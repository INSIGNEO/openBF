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

function _save_snapshot(obs::TUIObserver)
    fname = string(obs.sim_name, "_snapshot_", round(Int, obs.scalars.t), "s.txt")
    try
        open(fname, "w") do io
            println(io, "snapshot  t=$(obs.scalars.t)s  cycle=$(obs.scalars.passed_cycles)")
            for (i, site) in enumerate(obs.waveforms)
                p_vals = latest(site.P, _WAVEFORM_SAMPLES)
                isempty(p_vals) && continue
                println(io, "  $(site.label)  P_min=$(minimum(p_vals)) P_max=$(maximum(p_vals)) mmHg")
            end
        end
        push!(obs.log, (time(), "[Info] snapshot saved → $fname"))
    catch e
        push!(obs.log, (time(), "[Warn] snapshot failed: $(sprint(showerror, e))"))
    end
    nothing
end

function _process_commands!(obs::TUIObserver)
    while isready(obs.commands)
        cmd = take!(obs.commands)
        if cmd === :quit
            obs.should_stop[] = true
        elseif cmd === :pause
            obs.paused[] = !obs.paused[]
        elseif cmd === :tab
            n = length(_PANE_LABELS)
            obs.active_pane[] = (obs.active_pane[] % n) + 1
        elseif cmd === :snapshot
            _save_snapshot(obs)
        end
    end
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
            _process_commands!(obs)
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
    obs.paused[]      = false
    start_input!(obs)
    obs.render_task = if Threads.nthreads(:interactive) > 0
        Threads.@spawn :interactive _render_loop(obs)
    else
        Threads.@spawn :default _render_loop(obs)
    end
    nothing
end

function stop_render!(obs::TUIObserver)
    obs.should_stop[] = true
    stop_input!(obs)
    t = obs.render_task
    t !== nothing && timedwait(() -> istaskdone(t), 2.0; pollint = 0.05)
    print(stdout, _CURSOR_SHOW, _ALT_EXIT); flush(stdout)
    uninstall_log_sink!(obs)
    nothing
end
