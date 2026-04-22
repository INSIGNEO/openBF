using Term: Panel
using UnicodePlots: lineplot, lineplot!, BrailleCanvas

include("ringbuf.jl")
include("snapshot.jl")
include("observer.jl")
include("loghook.jl")
include("panes/waveforms.jl")
include("panes/header.jl")
include("panes/sidebar.jl")
include("panes/convergence.jl")
include("panes/log.jl")
include("renderer.jl")

function tui_should_run(cli_flag::Bool)
    cli_flag                              || return false
    isinteractive() || stdout isa Base.TTY || return false
    !haskey(ENV, "CI")                    || return false
    !haskey(ENV, "OPENBF_NO_TUI")         || return false
    return true
end
