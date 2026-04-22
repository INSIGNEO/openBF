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

function _render_frame(obs::TUIObserver, buf::IOBuffer)
    _update_rate!(obs)
    obs.frame_count += 1
    _, cols = _term_size()

    header  = _make_header(obs, cols)
    wave    = draw_waveforms(obs, cols)
    sidebar = _make_sidebar(obs, _BODY_HEIGHT)

    truncate(buf, 0); seek(buf, 0)
    print(buf, _SCREEN_HOME)
    if wave !== nothing
        print(buf, header / (wave * sidebar))
    else
        # no waveform data yet — show header + waiting message
        waiting = Panel(" Waiting for simulation data… (stride=$(obs.step_stride) steps)";
                        width = cols, padding = (0,0,0,0), style = "dim")
        print(buf, header / waiting)
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
                # show the error on screen so it's diagnosable
                msg = sprint(showerror, e)
                if msg != err_msg
                    err_msg = msg
                    truncate(buf, 0); seek(buf, 0)
                    print(buf, _SCREEN_HOME, _SCREEN_CLEAR)
                    print(buf, "TUI render error (solver still running):\n", msg, "\n")
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
    print(stdout, _ALT_ENTER, _CURSOR_HIDE); flush(stdout)
    obs.should_stop[] = false
    obs.render_task = Threads.@spawn :interactive _render_loop(obs)
    nothing
end

function stop_render!(obs::TUIObserver)
    obs.should_stop[] = true
    t = obs.render_task
    t !== nothing && timedwait(() -> istaskdone(t), 2.0; pollint = 0.05)
    print(stdout, _CURSOR_SHOW, _ALT_EXIT); flush(stdout)
    nothing
end
