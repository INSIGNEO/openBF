using Test
using openBF: TUIObserver, MonitorSpec, WaveformSite, ConvergenceHistory, Snapshot,
              LogRing, RingBuffer, record_step!, start_render!, stop_render!

struct MockVessel
    A::Vector{Float64}
    A0::Vector{Float64}
    beta::Vector{Float64}
    Pext::Float64
    Pout::Float64
    Q::Vector{Float64}
end

function _make_obs()
    monitored = [MonitorSpec(1, 2, "aorta")]
    waveforms = [WaveformSite("aorta")]
    now = time()
    TUIObserver(waveforms, ConvergenceHistory(), Snapshot(), LogRing(128),
                monitored, 1, 0, Threads.Atomic{Bool}(false), nothing,
                "test", 1.0, 5, now,
                Threads.Atomic{Int}(1), Threads.Atomic{Int}(1), 0, 0, now, 0.0,
                Channel{Symbol}(Inf), Threads.Atomic{Bool}(false), nothing)
end

@testset "render thread lifecycle" begin

    @testset "start and stop cleanly" begin
        obs     = _make_obs()
        vessels = [MockVessel([1e-4, 2e-4, 3e-4], [2e-4, 2e-4, 2e-4],
                              [1e6, 1e6, 1e6], 0.0, 0.0, [0.1, 0.2, 0.3])]
        for k in 1:50
            record_step!(obs, vessels, k * 1e-3, 1e-3)
        end

        start_render!(obs)
        @test obs.render_task !== nothing
        @test !istaskdone(obs.render_task)

        sleep(0.15)  # ~4 frames at 30 Hz
        stop_render!(obs)

        @test obs.should_stop[]
        @test istaskdone(obs.render_task)
        @test !istaskfailed(obs.render_task)
    end

    @testset "stop before any data is safe" begin
        obs = _make_obs()
        start_render!(obs)
        stop_render!(obs)
        @test istaskdone(obs.render_task)
    end

    @testset "double stop_render! is a no-op" begin
        obs = _make_obs()
        start_render!(obs)
        stop_render!(obs)
        @test_nowarn stop_render!(obs)
    end

end
