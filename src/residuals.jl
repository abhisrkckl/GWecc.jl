export residual_PQR,
    residual_A, residual_px, residual, residuals, residuals_px, waveform_and_residuals

function residual_PQR(ecc::Eccentricity, scu::SinCos)
    e = ecc.e
    su = scu.sinx
    cu = scu.cosx
    c2u = cu * cu - su * su

    P = ((e + (-2 + e * e) * cu) * su) / (1 - e * cu)
    Q = (sqrt(1 - e^2) * (e * cu - c2u)) / (1 - e * cu)
    R = e * su

    return P, Q, R
end

function residual_A(ecc::Eccentricity, phase::OrbitalPhase)
    s2ω = phase.sc2ω.sinx
    c2ω = phase.sc2ω.cosx

    P, Q, R = residual_PQR(ecc, phase.scu)

    A0 = R
    A1 = -P * s2ω + Q * c2ω
    A2 = P * c2ω + Q * s2ω

    return A0, A1, A2
end

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

    A0, A1, A2 = residual_A(e, phase)
    a0, a1, a2 = waveform_coeffs_c(proj)

    h0 = gw_amplitude(mass, n, e, dl)
    s0 = h0 / n.n

    sA = s0 * (a1 * A1 + a0 * A0)
    sB = s0 * (a2 * A2)

    s2ψ, c2ψ = proj.sc2ψ.sinx, proj.sc2ψ.cosx

    sp = c2ψ * sA - s2ψ * sB
    sx = s2ψ * sA + c2ψ * sB

    return sp, sx
end

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
        sp = sp - spE
        sx = sx - sxE
    end

    if PULSAR in terms
        dtp = dt + Δp
        spP, sxP = residual_px(mass, coeffs, l0p, proj, dl, true, dtp)
        sp = sp + spP
        sx = sx + sxP
    end

    return ap.Fp * sp + ap.Fx * sx
end

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
    dts = [redshifted_time_difference(tE, tref, z) for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)
    ap = AntennaPattern(psrpos, gwpos)

    psrterm = term == PULSAR
    delay = psrterm ? pulsar_term_delay(ap, dp, z) : Time(0.0)

    spxs = [residual_px(mass, coeffs, l0p, proj, dl, psrterm, dt + delay) for dt in dts]
    sps = first.(spxs) * (1 + z.z)
    sxs = last.(spxs) * (1 + z.z)

    return sps, sxs
end


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
    dts = [redshifted_time_difference(tE, tref, z) for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)
    ap = AntennaPattern(psrpos, gwpos)
    Δp = pulsar_term_delay(ap, dp, z)

    ss =
        [residual(mass, coeffs, l0p, proj, dl, ap, terms, Δp, dt) * (1 + z.z) for dt in dts]

    return ss
end

function waveform_and_residuals(
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

    hs = [waveform(mass, coeffs, l0p, proj, dl, ap, terms, Δp, dt) for dt in dts]
    ss =
        [residual(mass, coeffs, l0p, proj, dl, ap, terms, Δp, dt) for dt in dts] * (1 + z.z)

    return hs, ss
end
