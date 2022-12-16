export residual_from_components, residuals_from_components, residuals_components_𝒜

function residual_components_𝒜(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l_init::Angle,
    ap::AntennaPattern,
    dl::Distance,
    dt::Time,
)
    Fp, Fx = ap.Fp, ap.Fx

    n, e, l, g = evolve_orbit(coeffs, l_init, Angle(0.0), dt)
    phase = OrbitalPhase(mass, n, e, l, g)

    S = gw_amplitude(mass, n, e, dl) / n.n

    s2g, c2g = phase.sc2ω.sinx, phase.sc2ω.cosx
    P, Q, R = residual_PQR(e, phase.scu)

    𝒜1 = S * Fp * (Q * c2g - P * s2g)
    𝒜2 = S * Fp * (P * c2g + Q * s2g)
    𝒜3 = S * Fx * (Q * c2g - P * s2g)
    𝒜4 = S * Fx * (P * c2g + Q * s2g)
    𝒜5 = S * Fp * R
    𝒜6 = S * Fx * R

    return [𝒜1, 𝒜2, 𝒜3, 𝒜4, 𝒜5, 𝒜6]
end

function residual_component_coeffs_a(proj::ProjectionParams, term::Term)
    sgn = term == EARTH ? -1 : 1

    sc2γ0 = SinCos(Angle(term == EARTH ? proj.γ0 : proj.γp))
    s2γ0, c2γ0 = sc2γ0.sinx, sc2γ0.cosx
    s2ψ, c2ψ = proj.sc2ψ.sinx, proj.sc2ψ.cosx
    c0, c1, c2 = waveform_coeffs_c(proj)

    a1 = sgn * (-s2ψ * s2γ0 * c2 + c2ψ * c2γ0 * c1)
    a2 = sgn * (-s2ψ * c2γ0 * c2 - c2ψ * s2γ0 * c1)
    a3 = sgn * (c2ψ * s2γ0 * c2 + s2ψ * c2γ0 * c1)
    a4 = sgn * (c2ψ * c2γ0 * c2 - s2ψ * s2γ0 * c1)
    a5 = sgn * (c2ψ * c0)
    a6 = sgn * (s2ψ * c0)

    return [a1, a2, a3, a4, a5, a6]
end

function residual_from_components(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l_init::Angle,
    proj::ProjectionParams,
    dl::Distance,
    ap::AntennaPattern,
    term::Term,
    dt::Time,
)
    𝒜 = residual_components_𝒜(mass, coeffs, l_init, ap, dl, dt)
    a = residual_component_coeffs_a(proj, term)

    return dot(𝒜, a)
end

function residuals_from_components(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l_init::Angle,
    proj::ProjectionParams,
    dl::Distance,
    dp::Distance,
    psrpos::SkyLocation,
    gwpos::SkyLocation,
    z::Redshift,
    terms::Vector{Term},
    tref::Time,
    tEs::Vector{Time},
)
    dts = [redshifted_time_difference(tE, tref, z) for tE in tEs]

    res = zeros(length(tEs))

    coeffs = EvolvCoeffs(mass, n_init, e_init)
    ap = AntennaPattern(psrpos, gwpos)

    if EARTH in terms
        res =
            res + [
                residual_from_components(mass, coeffs, l_init, proj, dl, ap, EARTH, dt) *
                (1 + z.z) for dt in dts
            ]
    end

    if PULSAR in terms
        delay = pulsar_term_delay(ap, dp, z)
        res =
            res + [
                residual_from_components(
                    mass,
                    coeffs,
                    l_init,
                    proj,
                    dl,
                    ap,
                    PULSAR,
                    dt + delay,
                ) * (1 + z.z) for dt in dts
            ]
    end

    return res
end

function residuals_components_𝒜(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l_init::Angle,
    dl::Distance,
    dp::Distance,
    psrpos::SkyLocation,
    gwpos::SkyLocation,
    z::Redshift,
    term::Term,
    tref::Time,
    tEs::Vector{Time},
)
    dts = [redshifted_time_difference(tE, tref, z) for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)
    ap = AntennaPattern(psrpos, gwpos)

    psrterm = term == PULSAR
    delay = psrterm ? pulsar_term_delay(ap, dp, z) : Time(0.0)

    𝒜s = [residual_components_𝒜(mass, coeffs, l_init, ap, dl, dt + delay) for dt in dts]

    return [[𝒜[idx] for 𝒜 in 𝒜s] for idx = 1:6]
end
