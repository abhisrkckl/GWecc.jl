export residual_from_components, residuals_from_components, residuals_components_𝒜

"PTA signal component functions for fast likelihood."
function residual_components_𝒜(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l0p::InitPhaseParams,
    ap::AntennaPattern,
    term::Term,
    dt::Time,
)
    Fp, Fx = ap.Fp, ap.Fx

    l_init = term == EARTH ? l0p.l0 : l0p.lp

    n, e, l, g = evolve_orbit(coeffs, l_init, Angle(0.0), dt)
    phase = OrbitalPhase(mass, n, e, l, g)

    s2g, c2g = phase.sc2ω.sinx, phase.sc2ω.cosx
    P, Q, R = residual_PQR(e, phase.scu)

    c = gwres_amplitude_ratio(mass, coeffs.n_init, coeffs.e_init, n, e)

    𝒜1 = Fp * c * (Q * c2g - P * s2g)
    𝒜2 = Fp * c * (P * c2g + Q * s2g)
    𝒜3 = Fx * c * (Q * c2g - P * s2g)
    𝒜4 = Fx * c * (P * c2g + Q * s2g)
    𝒜5 = Fp * c * R
    𝒜6 = Fx * c * R

    return [𝒜1, 𝒜2, 𝒜3, 𝒜4, 𝒜5, 𝒜6]
end

"PTA signal coefficients as functions of projection parameters."
function residual_component_coeffs_a(proj::ProjectionParams, term::Term)
    sgn = term == EARTH ? 1 : -1

    sc2γ0 = SinCos(Angle(term == EARTH ? proj.γ0 : proj.γp))
    s2γ0, c2γ0 = sc2γ0.sinx, sc2γ0.cosx
    s2ψ, c2ψ = proj.sc2ψ.sinx, proj.sc2ψ.cosx
    c0, c1, c2 = waveform_coeffs_c(proj)

    S0 = proj.S0

    a1 = S0 * sgn * (-s2ψ * s2γ0 * c2 + c2ψ * c2γ0 * c1)
    a2 = S0 * sgn * (-s2ψ * c2γ0 * c2 - c2ψ * s2γ0 * c1)
    a3 = S0 * sgn * (c2ψ * s2γ0 * c2 + s2ψ * c2γ0 * c1)
    a4 = S0 * sgn * (c2ψ * c2γ0 * c2 - s2ψ * s2γ0 * c1)
    a5 = S0 * sgn * (c2ψ * c0)
    a6 = S0 * sgn * (s2ψ * c0)

    return [a1, a2, a3, a4, a5, a6]
end

"PTA signal computed from components"
function residual_from_components(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    ap::AntennaPattern,
    term::Term,
    dt::Time,
)
    𝒜 = residual_components_𝒜(mass, coeffs, l0p, ap, term, dt)
    a = residual_component_coeffs_a(proj, term)

    return dot(𝒜, a)
end

"PTA signal computed from components"
function residuals_from_components(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    dp::Distance,
    psrpos::SkyLocation,
    gwpos::SkyLocation,
    terms::Vector{Term},
    tref::Time,
    tEs::Vector{Time},
)
    dts = [tE - tref for tE in tEs]

    res = zeros(length(tEs))

    coeffs = EvolvCoeffs(mass, n_init, e_init)
    ap = AntennaPattern(psrpos, gwpos)

    if EARTH in terms
        res =
            res + [
                residual_from_components(mass, coeffs, l0p, proj, ap, EARTH, dt) for
                dt in dts
            ]
    end

    if PULSAR in terms
        delay = pulsar_term_delay(ap, dp)
        res =
            res + [
                residual_from_components(
                    mass,
                    coeffs,
                    l0p,
                    proj,
                    ap,
                    PULSAR,
                    dt + delay,
                ) for dt in dts
            ]
    end

    return res
end

"PTA signal components for the fast likelihood"
function residuals_components_𝒜(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l0p::InitPhaseParams,
    dp::Distance,
    psrpos::SkyLocation,
    gwpos::SkyLocation,
    term::Term,
    tref::Time,
    tEs::Vector{Time},
)
    dts = [tE - tref for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)
    ap = AntennaPattern(psrpos, gwpos)

    psrterm = term == PULSAR
    delay = psrterm ? pulsar_term_delay(ap, dp) : Time(0.0)

    𝒜s = [residual_components_𝒜(mass, coeffs, l0p, ap, term, dt + delay) for dt in dts]

    return [[𝒜[idx] for 𝒜 in 𝒜s] for idx = 1:6]
end
