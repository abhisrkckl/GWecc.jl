using GWecc
using Test
using FiniteDifferences

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
        @test ScaledTime(0.0).τ == 0.0
        @test ScaledTime(10.0).τ == 10.0
        @test_throws DomainError ScaledTime(Inf)

        for e in [eccmin.e / 4, eccmin.e, 0.1, 0.5, 0.9, eccmax.e, (eccmax.e + 1) / 2]
            @test e_from_τ_from_e(0.1) ≈ 0.1
        end
        # @test isapprox(e_from_τ_from_e(eccmin.e / 4), eccmin.e / 4)
        # @test isapprox(e_from_τ_from_e(eccmin.e), eccmin.e)
        # @test isapprox(e_from_τ_from_e(0.1), 0.1)
        # @test isapprox(e_from_τ_from_e(0.5), 0.5)
        # @test isapprox(e_from_τ_from_e(0.9), 0.9)
        # @test isapprox(e_from_τ_from_e(eccmax.e), eccmax.e)
        # @test isapprox(e_from_τ_from_e((eccmax.e + 1) / 2), (eccmax.e + 1) / 2)
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
        @test n0.n ≈ n_init.n atol = 1e-9
        @test e0.e ≈ e_init.e atol = 1e-9
        @test l0.θ ≈ l_init.θ atol = 1e-8
        @test γ0.θ ≈ γ_init.θ atol = 1e-9

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

    @testset "derivatives" begin
        numdiff = central_fdm(5, 1)

        mass = Mass(5000.0, 0.1)
        n = MeanMotion(1e-8)
        l = Angle(0.0)
        γ = Angle(0.0)

        for e in [0.1, 0.5, 0.9]
            ecc = Eccentricity(e)
            coeffs = EvolvCoeffs(mass, n, ecc)

            dτ_de_anl = derivative_dτ_de(ecc)
            dτ_de_num = numdiff(τ_from_e, e)
            @test dτ_de_anl ≈ dτ_de_num atol = 1e-9

            τ = τ_from_e(e)
            de_dτ_num = numdiff(e_from_τ, τ)
            de_dτ_anl = 1 / dτ_de_anl
            @test de_dτ_anl ≈ de_dτ_num atol = 1e-9

            κ = evolv_coeff_κ(mass, n, ecc)
            de_dt_anl1 = derivative_de_dt(mass, n, ecc)
            de_dt_anl2 = -κ * de_dτ_anl
            @test de_dt_anl1 ≈ de_dt_anl2

            dn_dt_anl = derivative_dn_dt(mass, n, ecc)
            dn_de_anl = dn_dt_anl / de_dt_anl1
            n_from_e_func = e1 -> n_from_e(coeffs, Eccentricity(e1)).n
            dn_de_num = numdiff(n_from_e_func, e)
            @test dn_de_anl ≈ dn_de_num

            dlbar_de_anl = derivative_dlbar_de(ecc)
            dlbar_de_num = numdiff(lbar_from_e, e)
            @test dlbar_de_anl ≈ dlbar_de_num atol = 1e-9

            α = evolv_coeff_α(mass, n, ecc)
            dlbar_dτ_anl = dlbar_de_anl * de_dτ_anl
            dl_dτ_anl = -α * dlbar_dτ_anl
            dτ_dt_anl = -κ
            dl_dt_anl2 = dl_dτ_anl * dτ_dt_anl
            dl_dt_anl1 = n.n
            @test dl_dt_anl2 ≈ dl_dt_anl1 atol = 1e-9

            dγbar_de_anl = derivative_dγbar_de(ecc)
            dγbar_de_num = numdiff(γbar_from_e, e)
            @test dγbar_de_anl ≈ dγbar_de_num atol = 1e-9

            β = evolv_coeff_β(mass, n, ecc)
            dγbar_dτ_anl = dγbar_de_anl * de_dτ_anl
            dγ_dτ_anl = -β * dγbar_dτ_anl
            dγ_dt_anl2 = dγ_dτ_anl * dτ_dt_anl
            k1 = advance_of_periastron_1PN(mass, n, ecc).k
            dγ_dt_anl1 = k1 * n.n
            @test dγ_dt_anl2 ≈ dγ_dt_anl1 atol = 1e-9

            dγbar2_de_anl = derivative_dγbar2_de(ecc, mass)
            γbar2_from_e_func = e1 -> γbar2_from_e(Eccentricity(e1), mass)
            dγbar2_de_num = numdiff(γbar2_from_e_func, e)
            @test dγbar2_de_anl ≈ dγbar2_de_num atol = 1e-9

            β2 = evolv_coeff_β2(mass, n, ecc)
            dγbar2_dτ_anl = dγbar2_de_anl * de_dτ_anl
            dγ2_dτ_anl = -β2 * dγbar2_dτ_anl
            dγ2_dt_anl2 = dγ2_dτ_anl * dτ_dt_anl
            k2 = advance_of_periastron_2PN(mass, n, ecc).k
            dγ2_dt_anl1 = k2 * n.n
            @test dγ2_dt_anl2 ≈ dγ2_dt_anl1 atol = 1e-9

            dγbar3_de_anl = derivative_dγbar3_de(ecc, mass)
            γbar3_from_e_func = e1 -> γbar3_from_e(Eccentricity(e1), mass)
            dγbar3_de_num = numdiff(γbar3_from_e_func, e)
            # This comparison is so imprecise because γbar3_from_e implements a
            # Pade approximant rather than the exact expression. This is OK because
            # this is a 3PN correction.
            @test dγbar3_de_anl ≈ dγbar3_de_num rtol = 1e-2

            β3 = evolv_coeff_β3(mass, n, ecc)
            dγbar3_dτ_anl = dγbar3_de_anl * de_dτ_anl
            dγ3_dτ_anl = -β2 * dγbar3_dτ_anl
            dγ3_dt_anl2 = dγ3_dτ_anl * dτ_dt_anl
            k3 = advance_of_periastron_3PN(mass, n, ecc).k
            dγ3_dt_anl1 = k3 * n.n
            @test dγ3_dt_anl2 ≈ dγ3_dt_anl1 atol = 1e-9
        end
    end

    @testset "mikkola" begin
        e = Eccentricity(0.1)
        l = Angle(0.0)
        @test l == mikkola(e, l)

        for e in [0.1, 0.5, 0.9]
            ecc = Eccentricity(e)
            for _l in -7.0:0.5:7.0
                l = Angle(_l)
                u = mikkola(ecc, l)
                @test l.θ ≈ kepler(ecc, u).θ
            end
        end
    end

    @testset "orbital phase" begin
        for θ in -7.0:0.5:7.0
            scθ = SinCos(Angle(θ))
            @test scθ.sinx^2 + scθ.cosx^2 ≈ 1.0
        end
    end
end
