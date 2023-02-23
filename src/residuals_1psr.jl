export ProjectionParams1psr,
    compute_1psr_params,
    waveform_amplitude_ratio,
    residual_amplitude_ratio,
    waveform_1psr,
    residual_1psr,
    residuals_1psr,
    waveform_and_residual_1psr,
    waveform_and_residuals_1psr

function ProjectionParams1psr(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    proj::ProjectionParams,
    dl::Distance,
    psrpos::SkyLocation,
    gwpos::SkyLocation,
)
    H0 = gw_amplitude(mass, n_init, e_init, dl)
    ap = AntennaPattern(psrpos, gwpos)
    α = sqrt(ap.Fp^2 + ap.Fx^2)
    n0 = n_init.n

    dψ = acos(ap.Fp / α) / 2
    ψ = proj.sc2ψ.x.θ / 2 + dψ

    ci = proj.cosι
    c2ψ = cos(2 * ψ)
    s2ψ = sin(2 * ψ)
    c2γ0 = cos(2 * proj.γ0)
    s2γ0 = sin(2 * proj.γ0)

    β1 = (1 + ci^2) * c2ψ * c2γ0 - 2 * ci * s2ψ * s2γ0
    β2 = -(1 + ci^2) * c2ψ * s2γ0 - 2 * ci * s2ψ * c2γ0
    β3 = (1 - ci^2) * c2ψ

    β = sqrt(β1^2 + β2^2 + β3^2)
    σ = acos(β3 / β)
    ρ = atan(β2 / β1)

    ζ0 = H0 * α * β / n0

    return ProjectionParams1psr(ζ0, σ, ρ)
end

function compute_1psr_params(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    dl::Distance,
    dp::Distance,
    psrpos::SkyLocation,
    gwpos::SkyLocation,
    tref::Time,
)
    l_init = l0p.l0

    ap = AntennaPattern(psrpos, gwpos)
    Δp = pulsar_term_delay(ap, dp)

    proj1 = ProjectionParams1psr(mass, n_init, e_init, proj, dl, psrpos, gwpos)

    return mass, n_init, e_init, l_init, proj1, Δp, tref
end

function waveform_amplitude_ratio(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    n_now::MeanMotion,
    e_now::Eccentricity,
)::Float64
    x_init = pn_param_x(mass, n_init, e_init).x
    x_now = pn_param_x(mass, n_now, e_now).x
    return x_now / x_init
end

function residual_amplitude_ratio(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    n_now::MeanMotion,
    e_now::Eccentricity,
)::Float64
    S_init = pn_param_x(mass, n_init, e_init).x / n_init.n
    S_now = pn_param_x(mass, n_now, e_now).x / n_now.n
    return S_now / S_init
end

function waveform_coeffs_b(proj::ProjectionParams1psr)
    sσ, cσ = proj.scσ.sinx, proj.scσ.cosx
    sρ, cρ = proj.scρ.sinx, proj.scρ.cosx
    return cσ, sσ * cρ, sσ * sρ
end

function waveform_1psr_term(
    mass::Mass,
    coeffs::EvolvCoeffs,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l_init::Angle,
    proj::ProjectionParams1psr,
    dt::Time,
)
    γ_init = Angle(0.0)

    n, e, l, γ = evolve_orbit(coeffs, l_init, γ_init, dt)
    phase = OrbitalPhase(mass, n, e, l, γ)

    hB0, hB1, hB2 = waveform_A(e, phase)
    b0, b1, b2 = waveform_coeffs_b(proj)

    C = waveform_amplitude_ratio(mass, n_init, e_init, n, e)

    Z0 = proj.ζ0 * n_init.n

    return Z0 * C * (b1 * hB1 + b2 * hB2 + b0 * hB0)
end

function residual_1psr_term(
    mass::Mass,
    coeffs::EvolvCoeffs,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l_init::Angle,
    proj::ProjectionParams1psr,
    dt::Time,
)
    γ_init = Angle(0.0)

    n, e, l, γ = evolve_orbit(coeffs, l_init, γ_init, dt)
    phase = OrbitalPhase(mass, n, e, l, γ)

    sB0, sB1, sB2 = residual_A(e, phase)
    b0, b1, b2 = waveform_coeffs_b(proj)

    c = residual_amplitude_ratio(mass, n_init, e_init, n, e)

    ζ0 = proj.ζ0

    return ζ0 * c * (b1 * sB1 + b2 * sB2 + b0 * sB0)
end

function waveform_and_residual_1psr_term(
    mass::Mass,
    coeffs::EvolvCoeffs,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l_init::Angle,
    proj::ProjectionParams1psr,
    dt::Time,
)
    γ_init = Angle(0.0)

    n, e, l, γ = evolve_orbit(coeffs, l_init, γ_init, dt)
    phase = OrbitalPhase(mass, n, e, l, γ)

    sB0, sB1, sB2 = residual_A(e, phase)
    hB0, hB1, hB2 = waveform_A(e, phase)
    b0, b1, b2 = waveform_coeffs_b(proj)

    C = waveform_amplitude_ratio(mass, n_init, e_init, n, e)
    c = residual_amplitude_ratio(mass, n_init, e_init, n, e)

    ζ0 = proj.ζ0
    Z0 = ζ0 * n_init.n

    h = Z0 * C * (b1 * hB1 + b2 * hB2 + b0 * hB0)
    s = ζ0 * c * (b1 * sB1 + b2 * sB2 + b0 * sB0)

    return h, s
end

"PTA signal for the single-pulsar case."
function waveform_1psr(
    mass::Mass,
    coeffs::EvolvCoeffs,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l_init::Angle,
    proj::ProjectionParams1psr,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    h = 0.0

    if EARTH in terms
        h += waveform_1psr_term(mass, coeffs, n_init, e_init, l_init, proj, dt)
    end

    if PULSAR in terms
        dtp = dt + Δp
        h += waveform_1psr_term(mass, coeffs, n_init, e_init, l_init, proj, dtp)
    end

    return h
end

"PTA signal for the single-pulsar case."
function residual_1psr(
    mass::Mass,
    coeffs::EvolvCoeffs,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l_init::Angle,
    proj::ProjectionParams1psr,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    R = 0.0

    if EARTH in terms
        R += residual_1psr_term(mass, coeffs, n_init, e_init, l_init, proj, dt)
    end

    if PULSAR in terms
        dtp = dt + Δp
        R += residual_1psr_term(mass, coeffs, n_init, e_init, l_init, proj, dtp)
    end

    return R
end

"PTA signal for the single-pulsar case."
function waveform_and_residual_1psr(
    mass::Mass,
    coeffs::EvolvCoeffs,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l_init::Angle,
    proj::ProjectionParams1psr,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    R = 0.0
    h = 0.0

    if EARTH in terms
        hE, RE =
            waveform_and_residual_1psr_term(mass, coeffs, n_init, e_init, l_init, proj, dt)
        h += hE
        R += RE
    end

    if PULSAR in terms
        dtp = dt + Δp
        hP, RP =
            waveform_and_residual_1psr_term(mass, coeffs, n_init, e_init, l_init, proj, dtp)
        h += hP
        R += RP
    end

    return h, R
end

"PTA signal for the single-pulsar case"
function waveform_1psr(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l_init::Angle,
    proj::ProjectionParams1psr,
    Δp::Time,
    terms::Vector{Term},
    tref::Time,
    tEs::Vector{Time},
)
    dts = [(tE - tref) for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)

    hs = [
        waveform_1psr(mass, coeffs, n_init, e_init, l_init, proj, terms, Δp, dt) for
        dt in dts
    ]

    return hs
end

"PTA signal for the single-pulsar case"
function residuals_1psr(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l_init::Angle,
    proj::ProjectionParams1psr,
    Δp::Time,
    terms::Vector{Term},
    tref::Time,
    tEs::Vector{Time},
)
    dts = [(tE - tref) for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)

    Rs = [
        residual_1psr(mass, coeffs, n_init, e_init, l_init, proj, terms, Δp, dt) for
        dt in dts
    ]

    return Rs
end

"PTA signal for the single-pulsar case"
function waveform_and_residuals_1psr(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l_init::Angle,
    proj::ProjectionParams1psr,
    Δp::Time,
    terms::Vector{Term},
    tref::Time,
    tEs::Vector{Time},
)
    dts = [(tE - tref) for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)

    hRs = [
        waveform_and_residual_1psr(
            mass,
            coeffs,
            n_init,
            e_init,
            l_init,
            proj,
            terms,
            Δp,
            dt,
        ) for dt in dts
    ]

    hs = first.(hRs)
    Rs = last.(hRs)

    return hs, Rs
end
