using GWecc
using Test
using FiniteDifferences
using LinearAlgebra
using UnPack
using NumericalIntegration
using Statistics


e_from_τ_from_e(ecc::Float64)::Float64 = e_from_τ(τ_from_e(Eccentricity(ecc))).e

@testset verbose = true "GWecc" begin
    @testset "parameters" begin
        @testset "mass" begin
            @test_throws DomainError Mass(4e-2, 0.2)
            @test_throws DomainError Mass(1.0, -1.0)
            @test_throws DomainError Mass(1.0, 0.0)
            @test_throws DomainError Mass(6e5, 0.2)
            @test_throws DomainError Mass(1.0, 0.26)
            @test Mass(1.0, 0.1).m == 1.0
            @test Mass(1.0, 0.25).η == 0.25
            @test Mass(1.0, 0.1).m > Mass(1.0, 0.1).Mch
        end

        @testset "time" begin
            @test (-Time(-1.0)).t == (-1.0 * Time(-1.0)).t
            @test (Time(-1.0) * -1).t == (-1.0 * Time(-1.0)).t
        end

        @testset "gammabar" begin
            @test_throws DomainError ScaledPeriastronAngle(0.01, 0.01, -43.0)
            @test_throws DomainError ScaledPeriastronAngle(0.01, 2.50, -44.0)
            @test_throws DomainError ScaledPeriastronAngle(0.10, 0.24, -44.0)
        end

        @testset "init phase params" begin
            @test_throws DomainError InitPhaseParams(-7.0, 1.0)
            @test_throws DomainError InitPhaseParams(1.0, 7.0)
        end

        @testset "projection params" begin
            @test_throws DomainError ProjectionParams(-1e-9, 1.0, 0.3, 1.0, 1.0)
            @test_throws DomainError ProjectionParams(1e-9, 4.0, 0.3, 1.0, 1.0)
            @test_throws DomainError ProjectionParams(1e-9, 1.0, 1.1, 1.0, 1.0)
            @test_throws DomainError ProjectionParams(1e-9, 1.0, 0.3, 4.0, 1.0)
            @test_throws DomainError ProjectionParams(1e-9, 1.0, 0.3, 1.0, 4.0)

            @test_throws DomainError ProjectionParams1psr(-1e-9, 0.1, 0.1)
            @test_throws DomainError ProjectionParams1psr(1e-9, -0.1, 0.1)
            @test_throws DomainError ProjectionParams1psr(1e-9, 0.1, 10.0)
        end

        @testset "sky location params" begin
            @test_throws DomainError SkyLocation(-1.0, 1.0)
            @test_throws DomainError SkyLocation(1.0, 4.1)
        end
    end

    @testset "orbital evlution" begin
        @testset "e and τ" begin
            @test_throws DomainError Eccentricity(0.0)
            @test Eccentricity(0.1).e == 0.1
            @test_throws DomainError Eccentricity(1.0)

            @test_throws DomainError ScaledTime(-1.0)
            @test ScaledTime(0.0).τ == 0.0
            @test ScaledTime(10.0).τ == 10.0
            @test_throws DomainError ScaledTime(Inf)

            for e in [eccmin.e / 4, eccmin.e, 0.1, 0.5, 0.9, eccmax.e, (eccmax.e + 1) / 2]
                @test e_from_τ_from_e(e) ≈ e
            end
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

    @testset "orbital phase" begin
        @testset "sincos" begin
            for θ = -7.0:0.5:7.0
                scθ = SinCos(Angle(θ))
                @test scθ.sinx^2 + scθ.cosx^2 ≈ 1.0
            end
        end

        @testset "mikkola" begin
            e = Eccentricity(0.1)
            l = Angle(0.0)
            @test l == mikkola(e, l)

            for e in [0.1, 0.5, 0.9]
                ecc = Eccentricity(e)
                for _l = -7.0:0.5:7.0
                    l = Angle(_l)
                    u = mikkola(ecc, l)
                    @test l.θ ≈ kepler(ecc, u).θ
                end
            end
        end

        @testset "true anomaly" begin
            for e in [0.1, 0.5, 0.9]
                ecc = Eccentricity(e)
                for _l in LinRange(-4 * π, 4 * π, 16)
                    l = Angle(_l)
                    u = mikkola(ecc, l)
                    scu = SinCos(u)
                    v_l = true_anomaly_diff(ecc, scu)
                    v = true_anomaly(ecc, scu)
                    @test v_l.θ ≈ v.θ - l.θ atol = 1e-9
                end
            end
        end

        mass = Mass(5000.0, 0.1)
        n = MeanMotion(1e-8)
        γ = Angle(3.5)
        γ_ = γ.θ

        for l_ in [0.0, 1.3, 2.3, 3.3, 4.3]
            l = Angle(l_)
            for e in [1e-9, 0.1, 0.5, 0.9]
                et = Eccentricity(e)
                orbital_phase = OrbitalPhase(mass, n, et, l, γ)
                φ = orbital_phase.sc2φ.x.θ / 2
                u = orbital_phase.scu.x

                if l == 0.0
                    @test φ ≈ γ_
                end

                if e == 1e-9
                    @test φ ≈ l_ + γ_
                end
            end
        end
    end

    @testset "antenna pattern" begin
        ra_p = 1.5
        dec_p = -0.8
        ra_gw = 0.5
        dec_gw = 0.75

        psrpos = SkyLocation(ra_p, dec_p)
        gwpos = SkyLocation(ra_gw, dec_gw)

        nhat = sky_direction_uvec(gwpos)
        @test norm(nhat) ≈ 1.0

        phat = sky_direction_uvec(psrpos)
        @test norm(phat) ≈ 1.0

        ep, ec = gw_polarization_tensors(gwpos)
        @test ep * ec + ec * ep ≈ zeros(3, 3) atol = 1e-9
        @test ep * nhat ≈ zeros(3) atol = 1e-9
        @test ec * nhat ≈ zeros(3) atol = 1e-9

        ap1 = AntennaPattern(psrpos, gwpos)
        ap2 = AntennaPattern(gwpos, psrpos)
        @test ap1.cosµ ≈ ap2.cosµ
        @test ap1.Fp^2 + ap1.Fx^2 ≈ ap2.Fp^2 + ap2.Fx^2

        @test_throws DomainError AntennaPattern(gwpos, gwpos)

        # @test AzimuthParam(0.5).α == 0.5
        # @test_throws DomainError AzimuthParam(1.1)
    end

    @testset "mismatch" begin
        a = [1, 2, 3]
        b = [2, -1, 0]
        K = [1 0 0; 0 1 0; 0 0 0]

        @test mismatch(a, a) ≈ 0.0
        @test mismatch(K, a, a) ≈ 0.0
        @test mismatch(a, -a) ≈ 2.0
        @test mismatch(K, a, -a) ≈ 2.0
        @test mismatch(a, b) ≈ 1.0
        @test mismatch(K, a, b) ≈ 1.0
    end

    @testset "waveform and residuals" begin
        mass = Mass(5000.0, 0.1)
        n_init = MeanMotion(1e-8)
        e_init = Eccentricity(0.1)
        l_init = Angle(0.1)
        γ_init = Angle(1.25)
        coeffs = EvolvCoeffs(mass, n_init, e_init)
        dl = Distance(1e15)

        h0 = gw_amplitude(mass, n_init, e_init, dl)
        @test h0 > 0

        n1 = MeanMotion(2e-8)
        e1 = Eccentricity(0.09)
        h1 = gw_amplitude(mass, n1, e1, dl)

        c01 = gwres_amplitude_ratio(mass, n_init, e_init, n1, e1)
        @test c01 ≈ (h1 / h0) / (n1.n / n_init.n)

        ra_p = 1.5
        dec_p = -0.8
        ra_gw = 0.5
        dec_gw = 0.75
        psrpos = SkyLocation(ra_p, dec_p)
        gwpos = SkyLocation(ra_gw, dec_gw)
        dp = Distance(1e13)
        ap = AntennaPattern(psrpos, gwpos)
        α = azimuth_param(ap)
        Δp = pulsar_term_delay(ap, dp)
        @test Δp.t < 0

        ψ = 1.1
        cosι = 0.52
        γ0 = γ_init.θ
        γp = γ0 + 0.2
        S0 = 1e-9
        proj = ProjectionParams(S0, ψ, cosι, γ0, γp)

        l0p = InitPhaseParams(l_init.θ, l_init.θ)

        dt = Time(10000.0)
        dtp = dt + Δp

        tref = Time(5 * 365.25 * 24 * 3600)
        tEs = time_range(Time(0.0), tref, 100)

        @testset "single point functions" begin

            sE = residual(mass, coeffs, l0p, proj, ap, [EARTH], Δp, dt)
            sP = residual(mass, coeffs, l0p, proj, ap, [PULSAR], Δp, dt)
            s = residual(mass, coeffs, l0p, proj, ap, [EARTH, PULSAR], Δp, dt)
            @test s ≈ sP + sE

            spE, sxE = residual_px(mass, coeffs, l0p, proj, false, dt)
            spP, sxP = residual_px(mass, coeffs, l0p, proj, true, dtp)
            @test sP ≈ -(ap.Fp * spP + ap.Fx * sxP)
            @test sE ≈ (ap.Fp * spE + ap.Fx * sxE)

            hE = waveform(mass, coeffs, l0p, proj, ap, [EARTH], Δp, dt)
            hP = waveform(mass, coeffs, l0p, proj, ap, [PULSAR], Δp, dt)
            h = waveform(mass, coeffs, l0p, proj, ap, [EARTH, PULSAR], Δp, dt)
            @test h ≈ hP + hE

            hpE, hxE = waveform_px(mass, coeffs, l0p, proj, false, dt)
            hpP, hxP = waveform_px(mass, coeffs, l0p, proj, true, dtp)
            @test hP ≈ -(ap.Fp * hpP + ap.Fx * hxP)
            @test hE ≈ (ap.Fp * hpE + ap.Fx * hxE)
        end

        @testset "terms and polarizations" begin
            rs = residuals(
                mass,
                n_init,
                e_init,
                l0p,
                proj,
                dp,
                psrpos,
                gwpos,
                [EARTH, PULSAR],
                tref,
                tEs,
            )
            rEs = residuals(
                mass,
                n_init,
                e_init,
                l0p,
                proj,
                dp,
                psrpos,
                gwpos,
                [EARTH],
                tref,
                tEs,
            )
            rPs = residuals(
                mass,
                n_init,
                e_init,
                l0p,
                proj,
                dp,
                psrpos,
                gwpos,
                [PULSAR],
                tref,
                tEs,
            )
            @test all(isfinite.(rs))
            @test all(isapprox.(rs, rEs + rPs))

            sps, sxs = residuals_px(
                mass,
                n_init,
                e_init,
                l0p,
                proj,
                dp,
                psrpos,
                gwpos,
                EARTH,
                tref,
                tEs,
            )
            @test all(isfinite.(sps)) && all(isfinite.(sxs))

            hs = waveform(
                mass,
                n_init,
                e_init,
                l0p,
                proj,
                dp,
                psrpos,
                gwpos,
                [EARTH, PULSAR],
                tref,
                tEs,
            )
            hEs = waveform(
                mass,
                n_init,
                e_init,
                l0p,
                proj,
                dp,
                psrpos,
                gwpos,
                [EARTH],
                tref,
                tEs,
            )
            hPs = waveform(
                mass,
                n_init,
                e_init,
                l0p,
                proj,
                dp,
                psrpos,
                gwpos,
                [PULSAR],
                tref,
                tEs,
            )
            @test all(isfinite.(hs))
            @test all(isapprox.(hs, hEs + hPs))

            rs1, hs1 = residuals_and_waveform(
                mass,
                n_init,
                e_init,
                l0p,
                proj,
                dp,
                psrpos,
                gwpos,
                [EARTH, PULSAR],
                tref,
                tEs,
            )
            @test all(isapprox.(hs1, hs))
            @test all(isapprox.(rs1, rs))
        end

        @testset "component functions" begin
            for e_init in Eccentricity.([0.1, 0.4, 0.8])
                sE = residuals(
                    mass,
                    n_init,
                    e_init,
                    l0p,
                    proj,
                    dp,
                    psrpos,
                    gwpos,
                    [EARTH],
                    tref,
                    tEs,
                )
                sP = residuals(
                    mass,
                    n_init,
                    e_init,
                    l0p,
                    proj,
                    dp,
                    psrpos,
                    gwpos,
                    [PULSAR],
                    tref,
                    tEs,
                )
                sEc = residuals_from_components(
                    mass,
                    n_init,
                    e_init,
                    l0p,
                    proj,
                    dp,
                    psrpos,
                    gwpos,
                    [EARTH],
                    tref,
                    tEs,
                )
                sPc = residuals_from_components(
                    mass,
                    n_init,
                    e_init,
                    l0p,
                    proj,
                    dp,
                    psrpos,
                    gwpos,
                    [PULSAR],
                    tref,
                    tEs,
                )
                @test all(isapprox.(sEc, sE, atol = 1e-8))
                @test all(isapprox.(sPc, sP, atol = 1e-8))

                𝒜s = residuals_components_𝒜(
                    mass,
                    n_init,
                    e_init,
                    l0p,
                    dp,
                    psrpos,
                    gwpos,
                    EARTH,
                    tref,
                    tEs,
                )
                @test all([all(isfinite.(𝒜)) for 𝒜 in 𝒜s])
            end
        end

        @testset "1psr functions" begin
            for e_init in Eccentricity.([0.1, 0.4, 0.8])
                dψ = polarization_angle_shift_1psr(ap)
                l0p0 = InitPhaseParams(l_init.θ)
                proj0 = ProjectionParams(S0, ψ, cosι, γ0, γ0)
                proj1 = ProjectionParams(S0 * α, ψ + dψ, cosι, γ0)
                proj1_ = ProjectionParams1psr(proj0, ap)
                Δp = pulsar_term_delay(ap, dp)
                ss = residuals(
                    mass,
                    n_init,
                    e_init,
                    l0p0,
                    proj0,
                    dp,
                    psrpos,
                    gwpos,
                    [EARTH, PULSAR],
                    tref,
                    tEs,
                )

                hs = waveform(
                    mass,
                    n_init,
                    e_init,
                    l0p0,
                    proj0,
                    dp,
                    psrpos,
                    gwpos,
                    [EARTH, PULSAR],
                    tref,
                    tEs,
                )

                hs_new = waveform_1psr(
                    mass,
                    n_init,
                    e_init,
                    l_init,
                    proj1_,
                    Δp,
                    [EARTH, PULSAR],
                    tref,
                    tEs,
                )
                @test all(isapprox.(hs, hs_new))

                ss_new = residuals_1psr(
                    mass,
                    n_init,
                    e_init,
                    l_init,
                    proj1_,
                    Δp,
                    [EARTH, PULSAR],
                    tref,
                    tEs,
                )
                @test all(isapprox.(ss, ss_new))
            end
        end

        @testset "h = ds/dt" begin
            year = 365.25 * 24 * 3600
            _ts = LinRange(0, 10 * year, 10000)
            tEs = Time.(_ts)
            tref = Time(maximum(_ts))

            for term in [EARTH, PULSAR]
                hps, hxs = waveform_px(
                    mass,
                    n_init,
                    e_init,
                    l0p,
                    proj,
                    dp,
                    psrpos,
                    gwpos,
                    term,
                    tref,
                    tEs,
                )
                rps, rxs = residuals_px(
                    mass,
                    n_init,
                    e_init,
                    l0p,
                    proj,
                    dp,
                    psrpos,
                    gwpos,
                    term,
                    tref,
                    tEs,
                )

                rps_n = cumul_integrate(_ts, hps)
                rxs_n = cumul_integrate(_ts, hxs)

                rps_n = rps_n .+ (mean(rps) - mean(rps_n))
                rxs_n = rxs_n .+ (mean(rxs) - mean(rxs_n))

                @test mismatch(rps, rps_n) < 1e-3
                @test mismatch(rxs, rxs_n) < 1e-3

                hs = waveform(
                    mass,
                    n_init,
                    e_init,
                    l0p,
                    proj,
                    dp,
                    psrpos,
                    gwpos,
                    [term],
                    tref,
                    tEs,
                )
                rs = residuals(
                    mass,
                    n_init,
                    e_init,
                    l0p,
                    proj,
                    dp,
                    psrpos,
                    gwpos,
                    [term],
                    tref,
                    tEs,
                )

                rs_n = cumul_integrate(_ts, hs)
                rs_n = rs_n .+ (mean(rs) - mean(rs_n))

                @test mismatch(rs, rs_n) < 1e-3
            end
        end

        @testset "spline functions" begin
            year = 365.25 * 24 * 3600
            _ts = LinRange(0, 10 * year, 10000)
            tEs = Time.(_ts)
            tref = Time(maximum(_ts))

            proj1 = ProjectionParams1psr(proj, ap)

            for term in [EARTH, PULSAR]
                rs = residuals(
                    mass,
                    n_init,
                    e_init,
                    l0p,
                    proj,
                    dp,
                    psrpos,
                    gwpos,
                    [term],
                    tref,
                    tEs,
                )
                rs_spl = residuals_spline(
                    mass,
                    n_init,
                    e_init,
                    l0p,
                    proj,
                    dp,
                    psrpos,
                    gwpos,
                    [term],
                    tref,
                    tEs,
                )
                @test mismatch(rs, rs_spl) < 1e-3

                Δp = pulsar_term_delay(ap, dp)
                rs = residuals_1psr(
                    mass,
                    n_init,
                    e_init,
                    l_init,
                    proj1,
                    Δp,
                    [term],
                    tref,
                    tEs,
                )
                rs_spl = residuals_1psr_spline(
                    mass,
                    n_init,
                    e_init,
                    l_init,
                    proj1,
                    Δp,
                    [term],
                    tref,
                    tEs,
                )
                @test mismatch(rs, rs_spl) < 1e-3
            end
        end
    end

    @testset "enterprise functions" begin
        toas = LinRange(0, 1000000, 100)
        pdist = 400.0 # kpc
        # alpha = azimuth_param(ap)
        psi = 1.1
        cos_inc = 0.5
        log10_M = 9.0
        eta = 0.25
        log10_F = -8.0
        e0 = 0.3
        gamma0 = 0.0
        gammap = 0.0
        sigma = 0.2
        rho = 0.15
        l0 = 0.0
        lp = 0.0
        tref = maximum(toas)
        log10_dl = -15.0
        deltap = 100.0
        for psrTerm in [true, false]
            res = eccentric_pta_signal_1psr(
                toas,
                sigma,
                rho,
                log10_M,
                eta,
                log10_F,
                e0,
                l0,
                tref,
                log10_dl,
                deltap,
                psrTerm,
            )
            @test all(isfinite.(res))
        end

        ra_p = 1.5
        dec_p = -0.8
        ra_gw = 0.5
        dec_gw = 0.75
        for psrTerm in [true, false]
            res = eccentric_pta_signal(
                toas,
                π / 2 - dec_p,
                ra_p,
                pdist,
                sin(dec_gw),
                ra_gw,
                psi,
                cos_inc,
                log10_M,
                eta,
                log10_F,
                e0,
                gamma0,
                gammap,
                l0,
                lp,
                tref,
                log10_dl,
                psrTerm,
            )
            @test all(isfinite.(res))
        end

        @test validate_params(log10_M, eta, log10_F, e0, tref, tref)
        @test !validate_params(10.0, eta, -7.0, 0.999, tref, tref)

        gwdist = 100.0
        log10_A = -9.0
        @test validate_params_target(log10_A, eta, log10_F, e0, gwdist, tref, tref)
        @test !validate_params_target(-5.0, eta, -7.0, 0.999, gwdist, tref, tref)
    end

end
