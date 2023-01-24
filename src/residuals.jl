export residual_PQR,
    residual_A,
    residual_px,
    residual,
    residual_1psr,
    residuals_px,
    residuals,
    residuals_and_waveform,
    residuals_1psr,
    residuals_and_waveform_1psr

"PTA signal component functions."
function residual_PQR(ecc::Eccentricity, scu::SinCos)
    e = ecc.e
    su = scu.sinx
    cu = scu.cosx
    c2u = cu * cu - su * su

    P = (sqrt(1 - e^2) * (c2u - e * cu)) / (1 - e * cu)
    Q = (((e^2 - 2) * cu + e) * su) / (1 - e * cu)
    R = e * su

    return P, Q, R
end

"PTA signal component functions."
function residual_A(ecc::Eccentricity, phase::OrbitalPhase)
    s2ω = phase.sc2ω.sinx
    c2ω = phase.sc2ω.cosx

    P, Q, R = residual_PQR(ecc, phase.scu)

    sA0 = R
    sA1 = -P * s2ω + Q * c2ω
    sA2 = P * c2ω + Q * s2ω

    return sA0, sA1, sA2
end

"+/x polarizations of the PTA signal."
function residual_px(
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

    sA0, sA1, sA2 = residual_A(e, phase)
    a0, a1, a2 = waveform_coeffs_c(proj)

    h0 = gw_amplitude(mass, n, e, dl)
    s0 = h0 / n.n

    sA = s0 * (a1 * sA1 + a0 * sA0)
    sB = s0 * (a2 * sA2)

    s2ψ, c2ψ = proj.sc2ψ.sinx, proj.sc2ψ.cosx

    sp = c2ψ * sA - s2ψ * sB
    sx = s2ψ * sA + c2ψ * sB

    return sp, sx
end

function residual_and_waveform_px(
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

    sA0, sA1, sA2 = residual_A(e, phase)
    hA0, hA1, hA2 = waveform_A(e, phase)
    a0, a1, a2 = waveform_coeffs_c(proj)

    h0 = gw_amplitude(mass, n, e, dl)
    s0 = h0 / n.n

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

"PTA signal."
function residual(
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
    sp = 0.0
    sx = 0.0

    if EARTH in terms
        spE, sxE = residual_px(mass, coeffs, l0p, proj, dl, false, dt)
        sp = sp + spE
        sx = sx + sxE
    end

    if PULSAR in terms
        dtp = dt + Δp
        spP, sxP = residual_px(mass, coeffs, l0p, proj, dl, true, dtp)
        sp = sp - spP
        sx = sx - sxP
    end

    return ap.Fp * sp + ap.Fx * sx
end

function residual_and_waveform(
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
    sp = 0.0
    sx = 0.0
    hp = 0.0
    hx = 0.0

    if EARTH in terms
        spE, sxE, hpE, hxE =
            residual_and_waveform_px(mass, coeffs, l0p, proj, dl, false, dt)
        sp = sp + spE
        sx = sx + sxE
        hp = hp + hpE
        hx = hx + hxE
    end

    if PULSAR in terms
        dtp = dt + Δp
        spP, sxP, hpP, hxP =
            residual_and_waveform_px(mass, coeffs, l0p, proj, dl, true, dtp)
        sp = sp - spP
        sx = sx - sxP
        hp = hp - hpP
        hx = hx - hxP
    end

    s = ap.Fp * sp + ap.Fx * sx
    h = ap.Fp * hp + ap.Fx * hx

    return s, h
end

"PTA signal for the single-pulsar case."
function residual_1psr(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    dl::Distance,
    α::AzimuthParam,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    sp = 0.0
    # sx = 0.0

    if EARTH in terms
        spE, sxE = residual_px(mass, coeffs, l0p, proj, dl, false, dt)
        sp = sp + spE
        # sx = sx + sxE
    end

    if PULSAR in terms
        dtp = dt + Δp
        spP, sxP = residual_px(mass, coeffs, l0p, proj, dl, true, dtp)
        sp = sp - spP
        # sx = sx - sxP
    end

    return α.α * sp
end

"PTA signal for the single-pulsar case."
function residual_and_waveform_1psr(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    dl::Distance,
    α::AzimuthParam,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    sp = 0.0
    hp = 0.0

    if EARTH in terms
        spE, sxE, hpE, hxE =
            residual_and_waveform_px(mass, coeffs, l0p, proj, dl, false, dt)
        sp = sp + spE
        hp = hp + hpE
    end

    if PULSAR in terms
        dtp = dt + Δp

        spP, sxP, hpP, hxP =
            residual_and_waveform_px(mass, coeffs, l0p, proj, dl, true, dtp)
        sp = sp - spP
        hp = hp - hpP
    end

    s = α.α * sp
    h = α.α * hp

    return s, h
end

"+/x polarizations of the PTA signal"
function residuals_px(
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
    term::Term,
    tref::Time,
    tEs::Vector{Time},
)
    dts = [unredshifted_time_difference(tE, tref, z) for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)
    ap = AntennaPattern(psrpos, gwpos)

    psrterm = term == PULSAR
    delay = psrterm ? pulsar_term_delay(ap, dp, z) : Time(0.0)

    spxs = [residual_px(mass, coeffs, l0p, proj, dl, psrterm, dt + delay) for dt in dts]
    sps = first.(spxs) * (1 + z.z)
    sxs = last.(spxs) * (1 + z.z)

    return sps, sxs
end

"PTA signal"
function residuals(
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
    dts = [unredshifted_time_difference(tE, tref, z) for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)
    ap = AntennaPattern(psrpos, gwpos)
    Δp = pulsar_term_delay(ap, dp, z)

    ss =
        [residual(mass, coeffs, l0p, proj, dl, ap, terms, Δp, dt) * (1 + z.z) for dt in dts]

    return ss
end

"PTA waveform and PTA signal"
function residuals_and_waveform(
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
    dts = [unredshifted_time_difference(tE, tref, z) for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)
    ap = AntennaPattern(psrpos, gwpos)
    Δp = pulsar_term_delay(ap, dp, z)

    shs = [
        residual_and_waveform(mass, coeffs, l0p, proj, dl, ap, terms, Δp, dt) for dt in dts
    ]

    ss = [sh[1] * (1 + z.z) for sh in shs]
    hs = [sh[2] for sh in shs]

    return ss, hs
end

"PTA signal for the single-pulsar case"
function residuals_1psr(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    dl::Distance,
    dp::Distance,
    α::AzimuthParam,
    z::Redshift,
    terms::Vector{Term},
    tref::Time,
    tEs::Vector{Time},
)
    dts = [unredshifted_time_difference(tE, tref, z) for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)
    Δp = pulsar_term_delay(α, dp, z)

    ss = [
        residual_1psr(mass, coeffs, l0p, proj, dl, α, terms, Δp, dt) * (1 + z.z) for
        dt in dts
    ]

    return ss
end

"PTA signal for the single-pulsar case"
function residuals_and_waveform_1psr(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    dl::Distance,
    dp::Distance,
    α::AzimuthParam,
    z::Redshift,
    terms::Vector{Term},
    tref::Time,
    tEs::Vector{Time},
)
    dts = [unredshifted_time_difference(tE, tref, z) for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)
    Δp = pulsar_term_delay(α, dp, z)

    shs = [
        residual_and_waveform_1psr(mass, coeffs, l0p, proj, dl, α, terms, Δp, dt) for
        dt in dts
    ]

    ss = first.(shs) * (1 + z.z)
    hs = last.(shs)

    return ss, hs
end
