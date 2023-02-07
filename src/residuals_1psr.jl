
function waveform_amplitude_ratio(
    mass::Mass,
    n_init::MeanMotion, 
    e_init::Eccentricity, 
    n_now::MeanMotion, 
    e_now::Eccentricity
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
    e_now::Eccentricity
)::Float64
    S_init = pn_param_x(mass, n_init, e_init).x / n_init.n
    S_now = pn_param_x(mass, n_now, e_now).x / n_now.n
    return S_now / S_init
end

function waveform_coeffs_b(proj::ProjectionParams1psr)
    sσ, cσ = proj.scσ.sinx, proj.scσ.cosx
    sρ, cρ = proj.scρ.sinx, proj.scρ.cosx
    return cσ, sσ*cρ, sσ*sρ
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
    b0, b1, b1 = waveform_coeffs_b(proj)

    C = waveform_amplitude_ratio(mass, n_init, e_init, n, e)
    
    Z0 = proj.ζ0 * n_init.n

    return Z0 * C * (b1*hB1 + b2*hB2 + b0*hB0)
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
    b0, b1, b1 = waveform_coeffs_b(proj)

    c = residual_amplitude_ratio(mass, n_init, e_init, n, e)
    
    ζ0 = proj.ζ0

    return ζ0 * c * (b1*sB1 + b2*sB2 + b0*sB0)
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
    b0, b1, b1 = waveform_coeffs_b(proj)

    C = waveform_amplitude_ratio(mass, n_init, e_init, n, e)
    c = residual_amplitude_ratio(mass, n_init, e_init, n, e)
    
    ζ0 = proj.ζ0
    Z0 = ζ0 * n_init.n

    h = Z0 * C * (b1*hB1 + b2*hB2 + b0*hB0)
    s = ζ0 * c * (b1*sB1 + b2*sB2 + b0*sB0)

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
        h += waveform_1psr_term(mass, coeffs, n_init, e_init, l_init, proj, dt)
        R += residual_1psr_term(mass, coeffs, n_init, e_init, l_init, proj, dt)
    end

    if PULSAR in terms
        dtp = dt + Δp
        h += waveform_1psr_term(mass, coeffs, n_init, e_init, l_init, proj, dtp)
        R += residual_1psr_term(mass, coeffs, n_init, e_init, l_init, proj, dtp)
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
    dts = [(tE-tref) for tE in tEs]

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
    dts = [(tE-tref) for tE in tEs]

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
    dts = [(tE-tref) for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)

    hRs = [
        waveform_and_residual_1psr(mass, coeffs, n_init, e_init, l_init, proj, terms, Δp, dt) for
        dt in dts
    ]

    hs = first.(hRs)
    Rs = last.(hRs)

    return hs, Rs
end