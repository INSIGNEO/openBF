# Golden-master regression oracle for all junction types.
#
# Models and the junction types they exercise:
#   ibif            — bifurcation (1 parent → 2 daughters)
#   adan56          — conjunction (1 parent → 1 daughter) + bifurcation
#   circle_of_willis — anastomosis (2 parents → 1 daughter)
#
# References serialised in test/ref_<model>.jls and committed to the repo.
# Any change that shifts a saved waveform by more than GOLDEN_TOL is a regression.

include("reference.jl")

const GOLDEN_TOL = 1e-12

const GOLDEN_MODELS = ("ibif", "adan56", "circle_of_willis")

@testset "Golden master" begin
    for model in GOLDEN_MODELS
        @testset "$model" begin
            ref  = load_reference(model)
            data, _ = capture_waveforms(model)
            for field in ("Q", "P", "A", "u")
                diffs = waveform_diff(ref, data; field=field)
                isempty(diffs) && continue
                @testset "$field" begin
                    for (label, d) in sort(collect(diffs))
                        @test d ≤ GOLDEN_TOL
                    end
                end
            end
        end
    end
end
