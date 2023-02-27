export gw_amplitude,
    gwres_amplitude_ratio, waveform_coeffs_c, waveform_px, waveform, waveform_1psr

"GW amplitude"
function gw_amplitude(
    mass::Mass,
    norb::MeanMotion,
    ecc::Eccentricity,
    dl::Distance,
)::Float64
    m, η = mass.m, mass.η
    dgw = dl.D
    x = pn_param_x(mass, norb, ecc).x
    return m * η * x / dgw
end

# function gw_amplitude_ratio(
#     mass::Mass,
#     n0::MeanMotion,
#     e0::Eccentricity,
#     n1::MeanMotion,
#     e1::Eccentricity,
# )
#     x0 = pn_param_x(mass, n0, e0).x
#     x1 = pn_param_x(mass, n1, e1).x
#     return x1 / x0
# end

function gwres_amplitude_ratio(
    mass::Mass,
    n0::MeanMotion,
    e0::Eccentricity,
    n1::MeanMotion,
    e1::Eccentricity,
)
    x0 = pn_param_x(mass, n0, e0).x
    x1 = pn_param_x(mass, n1, e1).x
    return (x1 / x0) / (n1.n / n0.n)
end

"Waveform coefficients that depend on the inclination"
function waveform_coeffs_c(proj::ProjectionParams)
    ci = proj.cosι
    return 1 - ci^2, 1 + ci^2, 2 * ci
end

"Waveform component functions"
function waveform_A(ecc::Eccentricity, phase::OrbitalPhase)
    e = ecc.e
    su = phase.scu.sinx
    cu = phase.scu.cosx
    s2φ = phase.sc2φ.sinx
    c2φ = phase.sc2φ.cosx

    χ = e * cu
    ξ = e * su

    P = (2 * e^2 - χ^2 + χ - 2) / (1 - χ)^2
    Q = (2 * sqrt(1 - e^2) * ξ) / (1 - χ)^2
    R = χ / (1 - χ)

    hA0 = R
    hA1 = -Q * s2φ + P * c2φ
    hA2 = Q * c2φ + P * s2φ

    return hA0, hA1, hA2
end

"+/x polarizations of the waveform."
function waveform_px(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    psrterm::Bool,
    dt::Time,
)
    γ_init = psrterm ? Angle(proj.γp) : Angle(proj.γ0)
    l_init = psrterm ? l0p.lp : l0p.l0

    n, e, l, γ = evolve_orbit(coeffs, l_init, γ_init, dt)
    phase = OrbitalPhase(mass, n, e, l, γ)

    hA0, hA1, hA2 = waveform_A(e, phase)
    a0, a1, a2 = waveform_coeffs_c(proj)

    c = gwres_amplitude_ratio(mass, coeffs.n_init, coeffs.e_init, n, e)
    s0 = proj.S0 * c
    h0 = s0 * n.n

    hA = h0 * (a1 * hA1 + a0 * hA0)
    hB = h0 * (a2 * hA2)

    s2ψ, c2ψ = proj.sc2ψ.sinx, proj.sc2ψ.cosx

    hp = c2ψ * hA - s2ψ * hB
    hx = s2ψ * hA + c2ψ * hB

    return hp, hx
end

"PTA waveform."
function waveform(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    ap::AntennaPattern,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    hp = 0.0
    hx = 0.0

    if EARTH in terms
        hpE, hxE = waveform_px(mass, coeffs, l0p, proj, false, dt)
        hp = hp + hpE
        hx = hx + hxE
    end

    if PULSAR in terms
        dtp = dt + Δp
        hpP, hxP = waveform_px(mass, coeffs, l0p, proj, true, dtp)
        hp = hp - hpP
        hx = hx - hxP
    end

    return ap.Fp * hp + ap.Fx * hx
end

"PTA waveform for single pulsar case."
function waveform_1psr(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l_init::Angle,
    proj::ProjectionParams,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    hp = 0.0
    # hx = 0.0

    l0p = InitPhaseParams(l_init.θ)

    if EARTH in terms
        hpE, hxE = waveform_px(mass, coeffs, l0p, proj, false, dt)
        hp = hp + hpE
        # hx = hx + hxE
    end

    if PULSAR in terms
        dtp = dt + Δp
        hpP, hxP = waveform_px(mass, coeffs, l0p, proj, true, dtp)
        hp = hp - hpP
        # hx = hx - hxP
    end

    return hp
end

"+/x polarizations of the PTA signal"
function waveform_px(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
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

    hpxs = [waveform_px(mass, coeffs, l0p, proj, psrterm, dt + delay) for dt in dts]
    hps = first.(hpxs)
    hxs = last.(hpxs)

    return hps, hxs
end

"PTA waveform"
function waveform(
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

    coeffs = EvolvCoeffs(mass, n_init, e_init)
    ap = AntennaPattern(psrpos, gwpos)
    Δp = pulsar_term_delay(ap, dp)

    ss = [waveform(mass, coeffs, l0p, proj, ap, terms, Δp, dt) for dt in dts]

    return ss
end

"PTA waveform for single pulsar case"
function waveform_1psr(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l_init::Angle,
    proj::ProjectionParams,
    Δp::Time,
    terms::Vector{Term},
    tref::Time,
    tEs::Vector{Time},
)
    dts = [tE - tref for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)

    ss = [waveform_1psr(mass, coeffs, l_init, proj, terms, Δp, dt) for dt in dts]

    return ss
end
