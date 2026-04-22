using Test
using openBF: TUIObserver, MonitorSpec, WaveformSite, ConvergenceHistory, Snapshot,
              LogRing, RingBuffer, record_step!, record_cycle!, latest

# Minimal duck-typed vessel for testing — only needs .P, .Q, .node3
struct MockVessel
    P::Vector{Float64}
    Q::Vector{Float64}
    node3::Int
end

function make_observer(n_sites::Int = 2; stride::Int = 10)
    monitored = [MonitorSpec(i, 2, "vessel_$i") for i in 1:n_sites]
    waveforms = [WaveformSite(spec.label) for spec in monitored]
    TUIObserver(waveforms, ConvergenceHistory(), Snapshot(), LogRing(128), monitored, stride, 0)
end

@testset "TUIObserver recording" begin

    @testset "waveform buffer sample count" begin
        stride = 10
        n_steps = 1000
        obs = make_observer(2; stride)
        vessels = [MockVessel([1.0, float(i), 3.0], [0.1, 0.2 * i, 0.3], 2) for i in 1:2]

        for k in 1:n_steps
            record_step!(obs, vessels, k * 1e-4, 1e-4)
        end

        expected = n_steps ÷ stride
        @test length(obs.waveforms[1].P) == expected
        @test length(obs.waveforms[2].Q) == expected
    end

    @testset "scalars updated" begin
        obs = make_observer(1; stride = 1)
        vessels = [MockVessel([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], 2)]
        record_step!(obs, vessels, 0.5, 1e-4)
        @test obs.scalars.t == 0.5
        @test obs.scalars.dt == 1e-4
    end

    @testset "convergence history updated" begin
        obs = make_observer(1; stride = 1)
        vessels = [MockVessel([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], 2)]
        record_cycle!(obs, 1, 0.05)
        record_cycle!(obs, 2, 0.02)
        vals = latest(obs.convergence.P, 10)
        @test length(vals) == 2
        @test vals[2] ≈ 0.02f0
    end

    @testset "stride decimation" begin
        stride = 7
        obs = make_observer(1; stride)
        vessels = [MockVessel([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], 2)]
        for k in 1:100
            record_step!(obs, vessels, k * 1e-4, 1e-4)
        end
        @test length(obs.waveforms[1].P) == 100 ÷ stride
    end

    @testset "zero allocations in hot path" begin
        obs = make_observer(4; stride = 1)
        vessels = [MockVessel(rand(10), rand(10), 5) for _ in 1:4]

        # warm up
        record_step!(obs, vessels, 0.0, 1e-4)
        obs.step_counter = 0  # reset so next call is a recording step

        allocs = @allocated record_step!(obs, vessels, 1e-4, 1e-4)
        @test allocs == 0
    end

end
