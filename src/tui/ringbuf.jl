mutable struct RingBuffer{T}
    data::Vector{T}
    capacity::Int
    write_idx::Threads.Atomic{Int}  # monotonically increasing count of total writes
end

function RingBuffer{T}(capacity::Int) where {T}
    RingBuffer{T}(Vector{T}(undef, capacity), capacity, Threads.Atomic{Int}(0))
end

function Base.push!(rb::RingBuffer{T}, x::T) where {T}
    k = Threads.atomic_add!(rb.write_idx, 1) + 1  # 1-based write count
    @inbounds rb.data[((k - 1) % rb.capacity) + 1] = x
    nothing
end

function latest(rb::RingBuffer{T}, n::Int) where {T}
    k = rb.write_idx[]
    n = min(n, k, rb.capacity)
    out = Vector{T}(undef, n)
    @inbounds for i in 1:n
        out[i] = rb.data[((k - n + i - 1) % rb.capacity) + 1]
    end
    out
end

Base.length(rb::RingBuffer) = min(rb.write_idx[], rb.capacity)
