export waveform_1psr,
    residual_1psr, residuals_1psr, residual_and_waveform_1psr, residuals_and_waveform_1psr

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

"PTA signal for the single-pulsar case."
function residual_1psr(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l_init::Angle,
    proj::ProjectionParams,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    sp = 0.0
    # sx = 0.0

    l0p = InitPhaseParams(l_init.θ)

    if EARTH in terms
        spE, sxE = residual_px(mass, coeffs, l0p, proj, false, dt)
        sp = sp + spE
        # sx = sx + sxE
    end

    if PULSAR in terms
        dtp = dt + Δp
        spP, sxP = residual_px(mass, coeffs, l0p, proj, true, dtp)
        sp = sp - spP
        # sx = sx - sxP
    end

    return sp
end

"PTA signal for the single-pulsar case"
function residuals_1psr(
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

    ss = [residual_1psr(mass, coeffs, l_init, proj, terms, Δp, dt) for dt in dts]

    return ss
end

"PTA signal for the single-pulsar case."
function residual_and_waveform_1psr(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l_init::Angle,
    proj::ProjectionParams,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    sp = 0.0
    hp = 0.0

    l0p = InitPhaseParams(l_init.θ)

    if EARTH in terms
        spE, sxE, hpE, hxE = residual_and_waveform_px(mass, coeffs, l0p, proj, false, dt)
        sp = sp + spE
        hp = hp + hpE
    end

    if PULSAR in terms
        dtp = dt + Δp

        spP, sxP, hpP, hxP = residual_and_waveform_px(mass, coeffs, l0p, proj, true, dtp)
        sp = sp - spP
        hp = hp - hpP
    end

    s = sp
    h = hp

    return s, h
end

"PTA signal for the single-pulsar case"
function residuals_and_waveform_1psr(
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

    shs = [
        residual_and_waveform_1psr(mass, coeffs, l_init, proj, terms, Δp, dt) for dt in dts
    ]

    ss = first.(shs)
    hs = last.(shs)

    return ss, hs
end
