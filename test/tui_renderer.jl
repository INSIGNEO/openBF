using Test
using openBF: TUIObserver, MonitorSpec, WaveformSite, ConvergenceHistory, Snapshot,
              LogRing, RingBuffer, record_step!, start_render!, stop_render!

struct MockVessel
    P::Vector{Float64}
    Q::Vector{Float64}
    node3::Int
end

@testset "render thread lifecycle" begin

    @testset "start and stop cleanly" begin
        monitored = [MonitorSpec(1, 2, "aorta")]
        waveforms = [WaveformSite("aorta")]
        obs = TUIObserver(waveforms, ConvergenceHistory(), Snapshot(), LogRing(128),
                          monitored, 1, 0, Threads.Atomic{Bool}(false), nothing)

        vessels = [MockVessel([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], 2)]
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
        monitored = [MonitorSpec(1, 2, "aorta")]
        waveforms = [WaveformSite("aorta")]
        obs = TUIObserver(waveforms, ConvergenceHistory(), Snapshot(), LogRing(128),
                          monitored, 1, 0, Threads.Atomic{Bool}(false), nothing)

        start_render!(obs)
        stop_render!(obs)
        @test istaskdone(obs.render_task)
    end

    @testset "double stop_render! is a no-op" begin
        monitored = [MonitorSpec(1, 2, "aorta")]
        waveforms = [WaveformSite("aorta")]
        obs = TUIObserver(waveforms, ConvergenceHistory(), Snapshot(), LogRing(128),
                          monitored, 1, 0, Threads.Atomic{Bool}(false), nothing)

        start_render!(obs)
        stop_render!(obs)
        @test_nowarn stop_render!(obs)
    end

end
