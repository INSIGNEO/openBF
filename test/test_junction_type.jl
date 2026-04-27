@testset "Junction type" begin
    @testset "construction — valid args" begin
        jc = Junction(7, [1, 2, 3], [:outlet, :inlet, :inlet])
        @test jc.id == 7
        @test jc.vessels == [1, 2, 3]
        @test jc.sides == [:outlet, :inlet, :inlet]
        @test jc.signs == [+1, -1, -1]
        @test jc.use_total_pressure == false
        @test length(jc.x)  == 6
        @test length(jc.F)  == 6
        @test size(jc.J)    == (6, 6)
        @test length(jc.dx) == 6
    end

    @testset "construction — k=2" begin
        jc = Junction(1, [4, 5], [:outlet, :inlet]; use_total_pressure=true)
        @test jc.use_total_pressure == true
        @test length(jc.x) == 4
        @test size(jc.J)   == (4, 4)
    end

    @testset "construction — k=4 (trifurcation)" begin
        jc = Junction(9, [1, 2, 3, 4], [:outlet, :inlet, :inlet, :inlet])
        @test jc.id == 9
        @test length(jc.vessels) == 4
        @test jc.signs == [+1, -1, -1, -1]
        @test jc.use_total_pressure == false
        @test length(jc.x)  == 8
        @test length(jc.F)  == 8
        @test size(jc.J)    == (8, 8)
        @test length(jc.dx) == 8
    end

    @testset "construction — invalid args" begin
        @test_throws ArgumentError Junction(1, [1, 2, 3], [:outlet, :inlet])
        @test_throws ArgumentError Junction(1, [1],       [:outlet])
        @test_throws ArgumentError Junction(1, [1, 2],    [:outlet, :bad])
    end

    @testset "construction — malformed sides" begin
        @test_throws AssertionError Junction(1, [1, 2], [:outlet, :outlet])
        @test_throws AssertionError Junction(1, [1, 2], [:inlet,  :inlet])
    end
end
