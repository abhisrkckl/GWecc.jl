export waveform, residual, residuals, residual_and_waveform, residuals_and_waveform

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

"PTA signal."
function residual(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    ap::AntennaPattern,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    sp = 0.0
    sx = 0.0

    if EARTH in terms
        spE, sxE = residual_px(mass, coeffs, l0p, proj, false, dt)
        sp = sp + spE
        sx = sx + sxE
    end

    if PULSAR in terms
        dtp = dt + Δp
        spP, sxP = residual_px(mass, coeffs, l0p, proj, true, dtp)
        sp = sp - spP
        sx = sx - sxP
    end

    return ap.Fp * sp + ap.Fx * sx
end

"PTA signal"
function residuals(
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

    ss = [residual(mass, coeffs, l0p, proj, ap, terms, Δp, dt) for dt in dts]

    return ss
end

function residual_and_waveform(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
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
        spE, sxE, hpE, hxE = residual_and_waveform_px(mass, coeffs, l0p, proj, false, dt)
        sp = sp + spE
        sx = sx + sxE
        hp = hp + hpE
        hx = hx + hxE
    end

    if PULSAR in terms
        dtp = dt + Δp
        spP, sxP, hpP, hxP = residual_and_waveform_px(mass, coeffs, l0p, proj, true, dtp)
        sp = sp - spP
        sx = sx - sxP
        hp = hp - hpP
        hx = hx - hxP
    end

    s = ap.Fp * sp + ap.Fx * sx
    h = ap.Fp * hp + ap.Fx * hx

    return s, h
end

"PTA waveform and PTA signal"
function residuals_and_waveform(
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

    shs = [residual_and_waveform(mass, coeffs, l0p, proj, ap, terms, Δp, dt) for dt in dts]

    ss = [sh[1] for sh in shs]
    hs = [sh[2] for sh in shs]

    return ss, hs
end
