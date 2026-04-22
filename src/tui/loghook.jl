using Logging: AbstractLogger, global_logger, LogLevel, Info
import Logging: handle_message, shouldlog, min_enabled_level, catch_exceptions

struct TUILogger <: AbstractLogger
    obs::TUIObserver
    prev::AbstractLogger
end

min_enabled_level(::TUILogger)              = Info
shouldlog(::TUILogger, args...)             = true
catch_exceptions(::TUILogger)               = false

function handle_message(l::TUILogger, level, msg, _module, group, id, file, line; kwargs...)
    push!(l.obs.log, (time(), string("[", level, "] ", msg)))
    nothing
end

function install_log_sink!(obs::TUIObserver)
    global_logger(TUILogger(obs, global_logger()))
    nothing
end
install_log_sink!(::Nothing) = nothing

function uninstall_log_sink!(obs::TUIObserver)
    lg = global_logger()
    lg isa TUILogger || return nothing
    global_logger(lg.prev)
    # flush the log ring to stdout now that the alt screen is gone
    entries = latest(obs.log, 128)
    if !isempty(entries)
        println("--- simulation log ($(length(entries)) entries) ---")
        for (t, msg) in entries
            @printf "  %6.0fs  %s\n" (t - obs.start_time) msg
        end
    end
    nothing
end
uninstall_log_sink!(::Nothing) = nothing
