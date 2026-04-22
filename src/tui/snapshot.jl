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

struct WaveformSite
    P::RingBuffer{Float32}
    Q::RingBuffer{Float32}
    label::String
end

WaveformSite(label::String, capacity::Int = 2048) =
    WaveformSite(RingBuffer{Float32}(capacity), RingBuffer{Float32}(capacity), label)

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
