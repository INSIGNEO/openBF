import REPL
using FileWatching: poll_fd

function _raw_mode!(enable::Bool)
    stdin isa Base.TTY || return
    try
        t = REPL.Terminals.TTYTerminal("", stdin, stdout, stderr)
        REPL.Terminals.raw!(t, enable)
    catch
    end
end

function _key_to_cmd(c::Char)
    c == 'q' || c == 'Q' || c == '\x03' ? :quit :   # \x03 = Ctrl-C in raw mode
    c == ' '                             ? :pause :
    c == '\t'                            ? :tab :
    c == 's' || c == 'S'                 ? :snapshot :
    nothing
end

function _input_loop(obs::TUIObserver)
    try
        _raw_mode!(true)
        while !obs.should_stop[]
            # 100 ms poll so the loop exits promptly when should_stop is set
            poll_fd(RawFD(0), 0.1; readable = true).readable || continue
            obs.should_stop[] && break
            c   = read(stdin, Char)
            cmd = _key_to_cmd(c)
            cmd !== nothing && put!(obs.commands, cmd)
        end
    catch e
        e isa InterruptException || rethrow()
    finally
        _raw_mode!(false)
    end
    nothing
end

function start_input!(obs::TUIObserver)
    stdin isa Base.TTY || return nothing   # skip on non-TTY stdin (pipes, CI)
    obs.input_task = if Threads.nthreads(:interactive) > 0
        Threads.@spawn :interactive _input_loop(obs)
    else
        Threads.@spawn :default _input_loop(obs)
    end
    nothing
end
start_input!(::Nothing) = nothing

function stop_input!(obs::TUIObserver)
    t = obs.input_task
    t !== nothing && timedwait(() -> istaskdone(t), 0.5; pollint = 0.05)
    nothing
end
stop_input!(::Nothing) = nothing
