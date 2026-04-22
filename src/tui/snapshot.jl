const _PA_TO_MMHG = Float32(1.0 / 133.322)

mutable struct Snapshot
    t::Float64
    dt::Float64
    cfl_vessel_id::Int
    cfl_margin::Float64
    cardiac_output::Float64
    map_pressure::Float64
    pulse_pressure::Float64
    gc_bytes::Int64
    passed_cycles::Int
end

Snapshot() = Snapshot(0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0)

# P is stored in mmHg (converted from Pa at record_cycle! time).
# curr_P / prev_P are full-cycle waveforms snapshotted from v.waveforms["P"][:, 4]
# once per cardiac cycle — no per-step pressure computation on the solver thread.
# Q is a live ring buffer of per-step samples for the sidebar stats.
mutable struct WaveformSite
    curr_P::Vector{Float32}
    prev_P::Vector{Float32}
    Q::RingBuffer{Float32}
    label::String
end

WaveformSite(label::String) =
    WaveformSite(Float32[], Float32[], RingBuffer{Float32}(2048), label)

struct ConvergenceHistory
    A::RingBuffer{Float32}
    Q::RingBuffer{Float32}
    P::RingBuffer{Float32}
end

ConvergenceHistory(capacity::Int = 512) =
    ConvergenceHistory(
        RingBuffer{Float32}(capacity),
        RingBuffer{Float32}(capacity),
        RingBuffer{Float32}(capacity),
    )

const LogRing = RingBuffer{Tuple{Float64,String}}
