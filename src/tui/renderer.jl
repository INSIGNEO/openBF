const _RENDER_HZ       = 30
const _RENDER_INTERVAL = 1.0 / _RENDER_HZ
const _CURSOR_HIDE     = "\e[?25l"
const _CURSOR_SHOW     = "\e[?25h"
const _ALT_ENTER       = "\e[?1049h"  # switch to alternate screen buffer
const _ALT_EXIT        = "\e[?1049l"  # restore normal screen buffer
const _SCREEN_HOME     = "\e[H"
const _SCREEN_CLEAR    = "\e[2J\e[H"

function _term_size()
    try
        displaysize(stdout)
    catch
        (24, 80)
    end
end

function _render_frame(obs::TUIObserver, buf::IOBuffer)
    truncate(buf, 0)
    seek(buf, 0)
    _, cols = _term_size()
    print(buf, _SCREEN_HOME)
    draw_waveforms(obs, buf, cols)
    write(stdout, take!(buf))
    flush(stdout)
    nothing
end

function _render_loop(obs::TUIObserver)
    buf = IOBuffer()
    # clear once on entry to the alternate screen
    print(stdout, _SCREEN_CLEAR)
    flush(stdout)
    try
        while !obs.should_stop[]
            sleep(_RENDER_INTERVAL)
            obs.should_stop[] && break
            try
                _render_frame(obs, buf)
            catch
            end
        end
    catch e
        e isa InterruptException || rethrow()
    end
    nothing
end

function start_render!(obs::TUIObserver)
    _prewarm_plots()
    print(stdout, _ALT_ENTER, _CURSOR_HIDE)
    flush(stdout)
    obs.should_stop[] = false
    obs.render_task = Threads.@spawn :interactive _render_loop(obs)
    nothing
end

function stop_render!(obs::TUIObserver)
    obs.should_stop[] = true
    t = obs.render_task
    if t !== nothing
        timedwait(() -> istaskdone(t), 2.0; pollint = 0.05)
    end
    print(stdout, _CURSOR_SHOW, _ALT_EXIT)
    flush(stdout)
    nothing
end
