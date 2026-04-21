using Test
using openBF: RingBuffer, latest

@testset "RingBuffer" begin

    @testset "basic push and latest" begin
        rb = RingBuffer{Int}(4)
        push!(rb, 10)
        push!(rb, 20)
        push!(rb, 30)
        @test latest(rb, 2) == [20, 30]
        @test latest(rb, 3) == [10, 20, 30]
    end

    @testset "write beyond capacity wraps" begin
        rb = RingBuffer{Int}(3)
        for i in 1:6
            push!(rb, i)
        end
        # only 3 slots; latest 3 should be 4,5,6
        @test latest(rb, 3) == [4, 5, 6]
        @test latest(rb, 5) == [4, 5, 6]  # clamps to capacity
    end

    @testset "latest n larger than written" begin
        rb = RingBuffer{Float32}(64)
        push!(rb, 1.0f0)
        push!(rb, 2.0f0)
        # requesting more than written returns what's there
        @test latest(rb, 10) == [1.0f0, 2.0f0]
    end

    @testset "length" begin
        rb = RingBuffer{Int}(8)
        @test length(rb) == 0
        push!(rb, 1); push!(rb, 2)
        @test length(rb) == 2
        for i in 3:12
            push!(rb, i)
        end
        @test length(rb) == 8  # capped at capacity
    end

    @testset "concurrent single-writer single-reader does not crash" begin
        rb = RingBuffer{Int}(256)
        n_writes = 10_000
        writer = Threads.@spawn begin
            for i in 1:n_writes
                push!(rb, i)
            end
        end
        # Reader just calls latest() while writer is running — must not crash or hang.
        reads = 0
        while rb.write_idx[] < n_writes
            latest(rb, 16)
            reads += 1
            yield()
        end
        wait(writer)
        # After writer done, final state must be the last `capacity` values written.
        vals = latest(rb, rb.capacity)
        @test length(vals) == rb.capacity
        @test maximum(vals) == n_writes
    end

end
