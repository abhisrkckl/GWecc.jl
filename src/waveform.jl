export gw_amplitude, waveform_coeffs_c, waveform_px, waveform

function gw_amplitude(mass::Mass, norb::MeanMotion, ecc::Eccentricity, dl::Distance)::Float64
    m, η = mass.m, mass.η
    dgw = dl.D
    x = pn_param_x(mass, norb, ecc).x
    return m * η * x / dgw
end

function waveform_coeffs_c(proj::ProjectionParams)
    ci = proj.cosι
    return 1 - ci^2, 1 + ci^2, 2 * ci
end

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

    A0 = R
    A1 = -Q * s2φ + P * c2φ
    A2 = Q * c2φ + P * s2φ

    return A0, A1, A2
end

function waveform_px(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    dl::Distance,
    psrterm::Bool,
    dt::Time,
)
    γ_init = psrterm ? Angle(proj.γp) : Angle(proj.γ0)
    l_init = psrterm ? l0p.lp : l0p.l0

    n, e, l, γ = evolve_orbit(coeffs, l_init, γ_init, dt)
    phase = OrbitalPhase(mass, n, e, l, γ)

    A0, A1, A2 = waveform_A(e, phase)
    a0, a1, a2 = waveform_coeffs_c(proj)

    h0 = gw_amplitude(mass, n, e, dl)

    hA = h0 * (a1 * A1 + a0 * A0)
    hB = h0 * (a2 * A2)

    s2ψ, c2ψ = proj.sc2ψ.sinx, proj.sc2ψ.cosx

    hp = c2ψ * hA - s2ψ * hB
    hx = s2ψ * hA + c2ψ * hB

    return hp, hx
end

function waveform(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    dl::Distance,
    ap::AntennaPattern,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    hp = 0.0
    hx = 0.0

    if EARTH in terms
        hpE, hxE = waveform_px(mass, coeffs, l0p, proj, dl, false, dt)
        hp = hp - hpE
        hx = hx - hxE
    end

    if PULSAR in terms
        dtp = dt + Δp
        hpP, hxP = waveform_px(mass, coeffs, l0p, proj, dl, true, dtp)
        hp = hp + hpP
        hx = hx + hxP
    end

    return ap.Fp * hp + ap.Fx * hx
end

function waveform(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l0p::InitPhaseParams,
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

    coeffs = EvolvCoeffs(mass, n_init, e_init)
    ap = AntennaPattern(psrpos, gwpos)
    Δp = pulsar_term_delay(ap, dp, z)

    ss = [waveform(mass, coeffs, l0p, proj, dl, ap, terms, Δp, dt) for dt in dts]

    return ss
end
