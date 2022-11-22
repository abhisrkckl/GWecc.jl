using GWecc
using Test

e_from_τ_from_e(ecc::Float64)::Float64 = e_from_τ(τ_from_e(Eccentricity(ecc))).e

@testset "GWecc" begin
    @testset "mass" begin
        @test_throws DomainError Mass(4e-2, 0.2)
        @test_throws DomainError Mass(1.0, -1.0)
        @test_throws DomainError Mass(1.0, 0.0)
        @test_throws DomainError Mass(6e4, 0.2)
        @test_throws DomainError Mass(1.0, 0.26)
        @test Mass(1.0, 0.1).m == 1.0
        @test Mass(1.0, 0.25).η == 0.25
    end

    @testset "e and τ" begin
        @test_throws DomainError Eccentricity(0.0)
        @test Eccentricity(0.1).e == 0.1
        @test_throws DomainError Eccentricity(1.0)

        @test_throws DomainError ScaledTime(-1.0)
        @test ScaledTime(10.0).τ == 10.0
        @test_throws DomainError ScaledTime(Inf)

        @test isapprox(e_from_τ_from_e(eccmin.e / 4), eccmin.e / 4)
        @test isapprox(e_from_τ_from_e(eccmin.e), eccmin.e)
        @test isapprox(e_from_τ_from_e(0.1), 0.1)
        @test isapprox(e_from_τ_from_e(0.5), 0.5)
        @test isapprox(e_from_τ_from_e(0.9), 0.9)
        @test isapprox(e_from_τ_from_e(eccmax.e), eccmax.e)
        @test isapprox(e_from_τ_from_e((eccmax.e + 1) / 2), (eccmax.e + 1) / 2)
    end

    @testset "n from e" begin
        @test_throws DomainError MeanMotion(6.2e-10)
        @test MeanMotion(1e-8).n == 1e-8
        @test_throws DomainError MeanMotion(6.4e-6)

        n0 = MeanMotion(1e-8)
        e0 = Eccentricity(0.1)
        e1 = Eccentricity(0.1)
        n1 = n_from_e(n0, e0, e1)
        @test isapprox(n1.n, n0.n)

        e1 = Eccentricity(0.2)
        n1 = n_from_e(n0, e0, e1)
        @test n1.n < n0.n

        e1 = Eccentricity(0.05)
        n1 = n_from_e(n0, e0, e1)
        @test n1.n > n0.n
    end

    @testset "coeffs" begin
        mass = Mass(1.0, 0.1)
        n_init = MeanMotion(1e-8)
        e_init = Eccentricity(0.1)
        coeffs = EvolvCoeffs(mass, n_init, e_init)
        @test coeffs.κ > 0 && isfinite(coeffs.κ)
        @test coeffs.α > 0 && isfinite(coeffs.α)
        @test coeffs.β > 0 && isfinite(coeffs.β)
        @test coeffs.β2 > 0 && isfinite(coeffs.β2)
        @test coeffs.β3 > 0 && isfinite(coeffs.β3)
    end

    @testset "lbar and γbar" begin
        mass = Mass(1.0, 0.1)
        n = MeanMotion(1e-8)
        e1 = Eccentricity(0.1)
        e2 = Eccentricity(0.2)

        lbar1 = lbar_from_e(e1)
        lbar2 = lbar_from_e(e2)
        @test lbar1.lbar < lbar2.lbar

        γbar1 = γbar_from_e(e1)
        γbar2 = γbar_from_e(e2)
        @test γbar1 < γbar2
    end
end
