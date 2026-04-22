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

# P: live rolling ring buffer in mmHg — per-step, for the progressive waveform display.
#    Pressure computation runs only when a TUIObserver is attached (Nothing dispatch is free).
# prev_P: previous cycle's full waveform in mmHg, snapshotted from v.waveforms["P"][:, 4]
#         once at cycle end — used as the comparison overlay.
# Q: live rolling ring buffer in m³/s — per-step, for sidebar stats.
mutable struct WaveformSite
    P::RingBuffer{Float32}
    prev_P::Vector{Float32}
    Q::RingBuffer{Float32}
    label::String
end

WaveformSite(label::String) =
    WaveformSite(RingBuffer{Float32}(2048), Float32[], RingBuffer{Float32}(2048), label)

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
