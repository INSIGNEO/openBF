include("observer.jl")
include("ringbuf.jl")
include("snapshot.jl")

function tui_should_run(cli_flag::Bool)
    cli_flag                              || return false
    isinteractive() || stdout isa Base.TTY || return false
    !haskey(ENV, "CI")                    || return false
    !haskey(ENV, "OPENBF_NO_TUI")         || return false
    return true
end
