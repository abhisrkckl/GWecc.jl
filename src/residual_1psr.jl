export waveform_1psr,
    residual_1psr, residuals_1psr, residual_and_waveform_1psr, residuals_and_waveform_1psr,
    residuals_1psr_new, waveform_1psr_new, residuals_and_waveform_1psr_new

# "PTA waveform for single pulsar case."
# function waveform_1psr(
#     mass::Mass,
#     coeffs::EvolvCoeffs,
#     l_init::Angle,
#     proj::ProjectionParams,
#     terms::Vector{Term},
#     Δp::Time,
#     dt::Time,
# )
#     hp = 0.0
#     # hx = 0.0

#     l0p = InitPhaseParams(l_init.θ)

#     if EARTH in terms
#         hpE, hxE = waveform_px(mass, coeffs, l0p, proj, false, dt)
#         hp = hp + hpE
#         # hx = hx + hxE
#     end

#     if PULSAR in terms
#         dtp = dt + Δp
#         hpP, hxP = waveform_px(mass, coeffs, l0p, proj, true, dtp)
#         hp = hp - hpP
#         # hx = hx - hxP
#     end

#     return hp
# end

# "PTA waveform for single pulsar case"
# function waveform_1psr(
#     mass::Mass,
#     n_init::MeanMotion,
#     e_init::Eccentricity,
#     l_init::Angle,
#     proj::ProjectionParams,
#     Δp::Time,
#     terms::Vector{Term},
#     tref::Time,
#     tEs::Vector{Time},
# )
#     dts = [tE - tref for tE in tEs]

#     coeffs = EvolvCoeffs(mass, n_init, e_init)

#     ss = [waveform_1psr(mass, coeffs, l_init, proj, terms, Δp, dt) for dt in dts]

#     return ss
# end

# "PTA signal for the single-pulsar case."
# function residual_1psr(
#     mass::Mass,
#     coeffs::EvolvCoeffs,
#     l_init::Angle,
#     proj::ProjectionParams,
#     terms::Vector{Term},
#     Δp::Time,
#     dt::Time,
# )
#     sp = 0.0
#     # sx = 0.0

#     l0p = InitPhaseParams(l_init.θ)

#     if EARTH in terms
#         spE, sxE = residual_px(mass, coeffs, l0p, proj, false, dt)
#         sp = sp + spE
#         # sx = sx + sxE
#     end

#     if PULSAR in terms
#         dtp = dt + Δp
#         spP, sxP = residual_px(mass, coeffs, l0p, proj, true, dtp)
#         sp = sp - spP
#         # sx = sx - sxP
#     end

#     return sp
# end

# "PTA signal for the single-pulsar case"
# function residuals_1psr(
#     mass::Mass,
#     n_init::MeanMotion,
#     e_init::Eccentricity,
#     l_init::Angle,
#     proj::ProjectionParams,
#     Δp::Time,
#     terms::Vector{Term},
#     tref::Time,
#     tEs::Vector{Time},
# )
#     dts = [tE - tref for tE in tEs]

#     coeffs = EvolvCoeffs(mass, n_init, e_init)

#     ss = [residual_1psr(mass, coeffs, l_init, proj, terms, Δp, dt) for dt in dts]

#     return ss
# end

# "PTA signal for the single-pulsar case."
# function residual_and_waveform_1psr(
#     mass::Mass,
#     coeffs::EvolvCoeffs,
#     l_init::Angle,
#     proj::ProjectionParams,
#     terms::Vector{Term},
#     Δp::Time,
#     dt::Time,
# )
#     sp = 0.0
#     hp = 0.0

#     l0p = InitPhaseParams(l_init.θ)

#     if EARTH in terms
#         spE, sxE, hpE, hxE = residual_and_waveform_px(mass, coeffs, l0p, proj, false, dt)
#         sp = sp + spE
#         hp = hp + hpE
#     end

#     if PULSAR in terms
#         dtp = dt + Δp

#         spP, sxP, hpP, hxP = residual_and_waveform_px(mass, coeffs, l0p, proj, true, dtp)
#         sp = sp - spP
#         hp = hp - hpP
#     end

#     s = sp
#     h = hp

#     return s, h
# end

# "PTA signal for the single-pulsar case"
# function residuals_and_waveform_1psr(
#     mass::Mass,
#     n_init::MeanMotion,
#     e_init::Eccentricity,
#     l_init::Angle,
#     proj::ProjectionParams,
#     Δp::Time,
#     terms::Vector{Term},
#     tref::Time,
#     tEs::Vector{Time},
# )
#     dts = [tE - tref for tE in tEs]

#     coeffs = EvolvCoeffs(mass, n_init, e_init)

#     shs = [
#         residual_and_waveform_1psr(mass, coeffs, l_init, proj, terms, Δp, dt) for dt in dts
#     ]

#     ss = first.(shs)
#     hs = last.(shs)

#     return ss, hs
# end

function residual_1psr_coeffs_β(proj::ProjectionParams)
    ci = proj.cosι
    c2ψ = proj.sc2ψ.cosx
    s2ψ = proj.sc2ψ.sinx
    c2γ0 = cos(2*proj.γ0)
    s2γ0 = sin(2*proj.γ0)

    β0 = (1-ci^2)*c2ψ
    β1 = (1+ci^2)*c2γ0*c2ψ - 2*ci*s2γ0*s2ψ
    β2 = -(1+ci^2)*s2γ0*c2ψ - 2*ci*c2γ0*s2ψ

    return β0, β1, β2
end

function ProjectionParams1psr(proj::ProjectionParams, ap::AntennaPattern)
    dψ = polarization_angle_shift_1psr(ap)
    α = AzimuthParam(ap).α
    ψ1 = proj.sc2ψ.x.θ/2 + dψ
    ci = proj.cosι
    γ0 = proj.γ0
    S0 = proj.S0
    proj1 = ProjectionParams(S0 * α, ψ1, ci, γ0)

    β0, β1, β2 = residual_1psr_coeffs_β(proj1)
    β = sqrt(β0^2 + β1^2 + β2^2)
    σ = atan(sqrt(β1^2 + β2^2), β0)
    ρ = atan(β2, β1)

    ζ0 = proj.S0 * α * β

    return ProjectionParams1psr(ζ0, σ, ρ)
end

function residual_1psr_coeffs_β(proj::ProjectionParams1psr)
    cσ = proj.scσ.cosx
    sσ = proj.scσ.sinx
    cρ = proj.scρ.cosx
    sρ = proj.scρ.sinx

    β0 = cσ
    β1 = sσ*cρ
    β2 = sσ*sρ

    return β0, β1, β2
end

function residual_1psr_term(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l_init::Angle,
    proj::ProjectionParams1psr,
    dt::Time,
)
    n, e, l, γ = evolve_orbit(coeffs, l_init, Angle(0.0), dt)
    phase = OrbitalPhase(mass, n, e, l, γ)

    sA0, sA1, sA2 = residual_A(e, phase)
    β0, β1, β2 = residual_1psr_coeffs_β(proj)

    ζ0 = proj.ζ0

    c = gwres_amplitude_ratio(mass, coeffs.n_init, coeffs.e_init, n, e)

    s = ζ0 * c * (β0*sA0 + β1*sA1 + β2*sA2)

    return s
end

function residual_1psr_new(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l_init::Angle,
    proj::ProjectionParams1psr,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    s = 0.0

    if EARTH in terms
        s += residual_1psr_term(mass, coeffs, l_init, proj, dt)
    end

    if PULSAR in terms
        dtp = dt + Δp
        s -= residual_1psr_term(mass, coeffs, l_init, proj, dtp)
    end

    return s
end

"PTA signal for the single-pulsar case"
function residuals_1psr_new(
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
    dts = [tE - tref for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)

    ss = [residual_1psr_new(mass, coeffs, l_init, proj, terms, Δp, dt) for dt in dts]

    return ss
end

function waveform_1psr_term(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l_init::Angle,
    proj::ProjectionParams1psr,
    dt::Time,
)
    n, e, l, γ = evolve_orbit(coeffs, l_init, Angle(0.0), dt)
    phase = OrbitalPhase(mass, n, e, l, γ)

    hA0, hA1, hA2 = waveform_A(e, phase)
    β0, β1, β2 = residual_1psr_coeffs_β(proj)

    ζ0 = proj.ζ0

    c = gwres_amplitude_ratio(mass, coeffs.n_init, coeffs.e_init, n, e)

    h = ζ0 * c * n.n * (β0*hA0 + β1*hA1 + β2*hA2)

    return h
end

function waveform_1psr_new(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l_init::Angle,
    proj::ProjectionParams1psr,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    h = 0.0

    if EARTH in terms
        h += waveform_1psr_term(mass, coeffs, l_init, proj, dt)
    end

    if PULSAR in terms
        dtp = dt + Δp
        h -= waveform_1psr_term(mass, coeffs, l_init, proj, dtp)
    end

    return h
end

"Waveform for the single-pulsar case"
function waveform_1psr_new(
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
    dts = [tE - tref for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)

    hs = [waveform_1psr_new(mass, coeffs, l_init, proj, terms, Δp, dt) for dt in dts]

    return hs
end

function residual_and_waveform_1psr_term(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l_init::Angle,
    proj::ProjectionParams1psr,
    dt::Time,
)
    n, e, l, γ = evolve_orbit(coeffs, l_init, Angle(0.0), dt)
    phase = OrbitalPhase(mass, n, e, l, γ)

    hA0, hA1, hA2 = waveform_A(e, phase)
    sA0, sA1, sA2 = residual_A(e, phase)
    β0, β1, β2 = residual_1psr_coeffs_β(proj)

    ζ0 = proj.ζ0

    c = gwres_amplitude_ratio(mass, coeffs.n_init, coeffs.e_init, n, e)

    h = ζ0 * c * n.n * (β0*hA0 + β1*hA1 + β2*hA2)
    s = ζ0 * c * (β0*sA0 + β1*sA1 + β2*sA2)

    return s, h
end

"PTA signal for the single-pulsar case."
function residual_and_waveform_1psr_new(
    mass::Mass,
    coeffs::EvolvCoeffs,
    l_init::Angle,
    proj::ProjectionParams1psr,
    terms::Vector{Term},
    Δp::Time,
    dt::Time,
)
    s = 0.0
    h = 0.0

    if EARTH in terms
        sE, hE = residual_and_waveform_1psr_term(mass, coeffs, l_init, proj, dt)
        s += sE
        h += hE
    end

    if PULSAR in terms
        dtp = dt + Δp
        sP, hP = residual_and_waveform_1psr_term(mass, coeffs, l_init, proj, dtp)
        s -= sP
        h -= hP
    end

    return s, h
end

"PTA signal for the single-pulsar case"
function residuals_and_waveform_1psr_new(
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
    dts = [tE - tref for tE in tEs]

    coeffs = EvolvCoeffs(mass, n_init, e_init)

    shs = [
        residual_and_waveform_1psr_new(mass, coeffs, l_init, proj, terms, Δp, dt) for dt in dts
    ]

    ss = first.(shs)
    hs = last.(shs)

    return ss, hs
end