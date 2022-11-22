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
        @test Mass(1.0, 0.1).m > Mass(1.0, 0.1).Mch
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

    @testset "evolve_orbit" begin
        mass = Mass(5000.0, 0.1)
        n_init = MeanMotion(1e-8)
        e_init = Eccentricity(0.1)
        l_init = Angle(0.1)
        γ_init = Angle(-1.25)
        coeffs = EvolvCoeffs(mass, n_init, e_init)

        delay0 = Time(0.0)
        n0, e0, l0, γ0 = evolve_orbit(coeffs, l_init, γ_init, delay0)

        @test n0.n ≈ n_init.n atol=1e-9
        @test e0.e ≈ e_init.e atol=1e-9
        @test_broken l0.θ ≈ l_init.θ atol=1e-9
        @test γ0.θ ≈ γ_init.θ atol=1e-9

        delay1 = Time(-10000.0)
        n1, e1, l1, γ1 = evolve_orbit(coeffs, l_init, γ_init, delay1)

        @test n1.n < n_init.n
        @test e1.e > e_init.e
        @test l1.θ < l_init.θ
        @test γ1.θ < γ_init.θ

        delay2 = Time(10000.0)
        n2, e2, l2, γ2 = evolve_orbit(coeffs, l_init, γ_init, delay2)

        @test n2.n > n_init.n
        @test e2.e < e_init.e
        @test l2.θ > l_init.θ
        @test γ2.θ > γ_init.θ
    end
end
