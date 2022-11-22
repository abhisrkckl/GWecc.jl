using GWecc
using Test

e_from_τ_from_e(ecc::Float64) :: Float64 = e_from_τ(τ_from_e(Eccentricity(ecc))).e

@testset "GWecc" begin
    @testset "parameters" begin
        @test_throws DomainError Eccentricity(0.0)
        @test Eccentricity(0.1).e == 0.1
        @test_throws DomainError Eccentricity(1.0)
        
        @test_throws DomainError ScaledTime(-1.0)
        @test ScaledTime(10.0).τ == 10.0
        @test_throws DomainError ScaledTime(Inf)
        
        @test_throws DomainError Mass(4e-2)
        @test Mass(1.0).m == 1.0
        @test_throws DomainError Mass(6e4)

        @test_throws DomainError MeanMotion(6.2e-10)
        @test MeanMotion(1e-8).n == 1e-8
        @test_throws DomainError MeanMotion(6.4e-6)
    end

    @testset "e and τ" begin    
        @test isapprox(e_from_τ_from_e(eccmin.e/4), eccmin.e/4)
        @test isapprox(e_from_τ_from_e(eccmin.e), eccmin.e)
        @test isapprox(e_from_τ_from_e(0.1), 0.1)
        @test isapprox(e_from_τ_from_e(0.5), 0.5)
        @test isapprox(e_from_τ_from_e(0.9), 0.9)
        @test isapprox(e_from_τ_from_e(eccmax.e), eccmax.e)
        @test isapprox(e_from_τ_from_e((eccmax.e+1)/2), (eccmax.e+1)/2)
    end


    
end