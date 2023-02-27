export waveform_px, residual_px, residuals_px

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

"+/x polarizations of the PTA signal."
function residual_px(
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

    sA0, sA1, sA2 = residual_A(e, phase)
    a0, a1, a2 = waveform_coeffs_c(proj)

    # h0 = gw_amplitude(mass, n, e, dl)
    c = gwres_amplitude_ratio(mass, coeffs.n_init, coeffs.e_init, n, e)
    s0 = proj.S0 * c

    sA = s0 * (a1 * sA1 + a0 * sA0)
    sB = s0 * (a2 * sA2)

    s2ψ, c2ψ = proj.sc2ψ.sinx, proj.sc2ψ.cosx

    sp = c2ψ * sA - s2ψ * sB
    sx = s2ψ * sA + c2ψ * sB

    return sp, sx
end

"+/x polarizations of the PTA signal"
function residuals_px(
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

    spxs = [residual_px(mass, coeffs, l0p, proj, psrterm, dt + delay) for dt in dts]
    sps = first.(spxs)
    sxs = last.(spxs)

    return sps, sxs
end

function residual_and_waveform_px(
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

    sA0, sA1, sA2 = residual_A(e, phase)
    hA0, hA1, hA2 = waveform_A(e, phase)
    a0, a1, a2 = waveform_coeffs_c(proj)

    c = gwres_amplitude_ratio(mass, coeffs.n_init, coeffs.e_init, n, e)
    s0 = proj.S0 * c
    h0 = s0 * n.n

    sA = s0 * (a1 * sA1 + a0 * sA0)
    sB = s0 * (a2 * sA2)

    hA = h0 * (a1 * hA1 + a0 * hA0)
    hB = h0 * (a2 * hA2)

    s2ψ, c2ψ = proj.sc2ψ.sinx, proj.sc2ψ.cosx

    sp = c2ψ * sA - s2ψ * sB
    sx = s2ψ * sA + c2ψ * sB

    hp = c2ψ * hA - s2ψ * hB
    hx = s2ψ * hA + c2ψ * hB

    return sp, sx, hp, hx
end

