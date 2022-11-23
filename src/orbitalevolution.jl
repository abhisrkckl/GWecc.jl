export e_from_τ, τ_from_e, eccmin, eccmax, taumin, taumax
export evolv_coeff_κ, evolv_coeff_α, evolv_coeff_β, evolv_coeff_β2, evolv_coeff_β3
export n_from_e, lbar_from_e, γbar_from_e, γbar2_from_e, γbar3_from_e
export EvolvCoeffs, evolve_orbit
export derivative_dτ_de, derivative_de_dt, derivative_dn_dt
export derivative_dlbar_de, derivative_dγbar_de, derivative_dγbar2_de, derivative_dγbar3_de
export advance_of_periastron_1PN, advance_of_periastron_2PN, advance_of_periastron_3PN

import DataInterpolations
import JLD
using HypergeometricFunctions

include("parameters.jl")

function read_precomputed_tau_e(datafile::String)
    data = JLD.load(datafile)
    return data["taus"], data["es"]
end

taus, es = read_precomputed_tau_e("/home/susobhan/Work/GWecc.jl/data/tau_e.jld")
tau_from_e_spline = DataInterpolations.CubicSpline(taus, es)
e_from_tau_spline = DataInterpolations.CubicSpline(es, taus)

eccmin = Eccentricity(minimum(es))
eccmax = Eccentricity(maximum(es))
taumin = ScaledTime(minimum(taus))
taumax = ScaledTime(maximum(taus))

a::Float64 = 2 * sqrt(2) / 5 / 5^(63 / 2299) / 17^(1181 / 2299);
b::Float64 = 2 / sqrt(1 - eccmax.e) - a * taumax.τ;

function e_from_τ(tau::ScaledTime)::Eccentricity
    τ = tau.τ
    τmin = taumin.τ
    τmax = taumax.τ

    if τ < τmin
        coeff::Float64 = 2^(559 / 726) * 3^(19 / 48) / 19^(145 / 242)
        e = coeff * τ^(19 / 48)
    elseif τ > τmax
        e = 1 - 4 / (a * τ + b)^2
    else
        e = e_from_tau_spline(τ)
    end

    return Eccentricity(e)
end
e_from_τ(τ::Float64)::Float64 = e_from_τ(ScaledTime(τ)).e

function τ_from_e(ecc::Eccentricity)::ScaledTime
    e = ecc.e
    emin = eccmin.e
    emax = eccmax.e

    if e < emin
        coeff = 19 * 19^(1181 / 2299) / (6 * 2^(2173 / 2299))
        τ = coeff * e^(48 / 19)
    elseif e > emax
        τ = (2 / sqrt(1 - e) - b) / a
    else
        τ = tau_from_e_spline(e)
    end

    return ScaledTime(τ)
end
τ_from_e(e::Float64)::Float64 = τ_from_e(Eccentricity(e)).τ

function evolv_coeff_κ(mass::Mass, n_init::MeanMotion, e_init::Eccentricity)::Float64
    Mch = mass.Mch
    n0 = n_init.n
    e0 = e_init.e

    return (1/15) * (Mch*n0)^(5 / 3) * n0 * e0^(48 / 19) * (304 + 121 * e0^2)^(3480 / 2299) /
           (1 - e0^2)^4
end

function evolv_coeff_α(mass::Mass, n_init::MeanMotion, e_init::Eccentricity)::Float64
    Mch = mass.Mch
    n0 = n_init.n
    e0 = e_init.e

    return 15 / Mch^(5 / 3) / n0^(5 / 3) / e0^(30 / 19) / (304 + 121 * e0^2)^(2175 / 2299) *
           (1 - e0^2)^2.5
end

function evolv_coeff_β(mass::Mass, n_init::MeanMotion, e_init::Eccentricity)::Float64
    Mch = mass.Mch
    M = mass.m
    n0 = n_init.n
    e0 = e_init.e

    return 45 / Mch^(5 / 3) * M^(2 / 3) / n0 / e0^(18 / 19) /
           (304 + 121 * e0^2)^(1305 / 2299) * (1 - e0^2)^1.5
end

function evolv_coeff_β2(mass::Mass, n_init::MeanMotion, e_init::Eccentricity)::Float64
    Mch = mass.Mch
    M = mass.m
    n0 = n_init.n
    e0 = e_init.e

    return 15 / (4 * Mch^(5 / 3)) * M^(4 / 3) / cbrt(n0) / e0^(6 / 19) /
           (304 + 121 * e0^2)^(435 / 2299) * sqrt(1 - e0^2)
end

function evolv_coeff_β3(mass::Mass, n_init::MeanMotion, e_init::Eccentricity)::Float64
    Mch = mass.Mch
    M = mass.m
    n0 = n_init.n
    e0 = e_init.e

    return 15 / (128 * Mch^(5 / 3)) *
           (M * M) *
           cbrt(n0) *
           e0^(6 / 19) *
           (304 + 121 * e0^2)^(435 / 2299) / sqrt(1 - e0^2)
end

function lbar_from_e(ecc::Eccentricity)::ScaledMeanAnomaly
    e = ecc.e

    coeff::Float64 = 19^(2175 / 2299) / 30 / 2^(496 / 2299)
    lbar = coeff * e * e^(11 / 19) * _₂F₁(124 / 2299, 15 / 19, 34 / 19, -121 * e^2 / 304)
    return ScaledMeanAnomaly(lbar)
end
lbar_from_e(e::Float64)::Float64 = lbar_from_e(Eccentricity(e)).lbar

function γbar_from_e(ecc::Eccentricity)::Float64
    e = ecc.e
    coeff::Float64 = 19^(1305 / 2299) / 36 / 2^(1677 / 2299)
    γbar = coeff * e^(18 / 19) * _₂F₁(994 / 2299, 9 / 19, 28 / 19, -121 * e^2 / 304)
    return γbar
end
γbar_from_e(e::Float64)::Float64 = γbar_from_e(Eccentricity(e))

function γbar2_from_e(ecc::Eccentricity, mass::Mass)::Float64
    e = ecc.e
    η = mass.η
    coeff::Float64 = 3 * 2^(1740 / 2299) * 19^(435 / 2299)
    γbar2 =
        e^(6 / 19) / 336 * (
            4 * (51 - 26 * η) * (304 + 121 * e^2)^(435 / 2299) +
            coeff * (23 + 2 * η) * _₂F₁(3 / 19, 1864 / 2299, 22 / 19, -121 * e^2 / 304)
        )
    return γbar2
end

function γbar3_from_e(ecc::Eccentricity, mass::Mass)::Float64
    #=
    γbar3(e) = (γbar3_0 + η*γbar3_1(e) + (η^2)*γbar3_2(e)) / e^(6/19)
     
    The functions gbar30 and gbar31 contain Appell F1 functions that are
    hard to compute. So I am using Pade approximants here. Since this is a
    3PN correction the error due to this approximation should be negligible.
    This probably won't work if e == 0.
    =#
    e = ecc.e
    η = mass.η
    γbar3_0::Float64 =
        -71.19148913450734 +
        (
            14.212810668086124 * (e * e) - 2.7809293422212327 * (e^4) -
            0.5247248712090506 * (e^6)
        ) / (1.0 + 0.021216346439953203 * (e * e) - 0.08582538305359699 * (e^4))
    γbar3_1::Float64 =
        75.17536852020817 +
        (
            -11.720605504532221 * (e * e) +
            0.9937219102426796 * (e^4) +
            0.3416818503954963 * (e^6)
        ) / (1.0 + 0.11825867242669967 * (e * e) - 0.05881762636038525 * (e^4))
    γbar3_2::Float64 =
        -3.1640661837558817 +
        (
            3.109257298381722 * (e * e) +
            1.0709804529990066 * (e^4) +
            0.05588268809419086 * (e^6)
        ) / (1.0 + 0.46114315553206664 * (e * e) + 0.041157522000398405 * (e^4))
    γbar3 = (γbar3_0 + η * γbar3_1 + η * η * γbar3_2) / e^(6 / 19)
    return γbar3
end

function γbar123_from_e(ecc::Eccentricity, mass::Mass)::ScaledPeriastronAngle
    γbar1 = γbar_from_e(ecc)
    γbar2 = γbar2_from_e(ecc, mass)
    γbar3 = γbar3_from_e(ecc, mass)
    return ScaledPeriastronAngle(γbar1, γbar2, γbar3)
end

struct EvolvCoeffs
    mass::Mass

    n_init::MeanMotion
    e_init::Eccentricity

    κ::Float64
    τ0::ScaledTime

    α::Float64
    lbar0::ScaledMeanAnomaly

    β::Float64
    β2::Float64
    β3::Float64
    γbar0_123::ScaledPeriastronAngle

    function EvolvCoeffs(mass::Mass, n_init::MeanMotion, e_init::Eccentricity)
        new(
            mass,
            n_init,
            e_init,
            evolv_coeff_κ(mass, n_init, e_init),
            τ_from_e(e_init),
            evolv_coeff_α(mass, n_init, e_init),
            lbar_from_e(e_init),
            evolv_coeff_β(mass, n_init, e_init),
            evolv_coeff_β2(mass, n_init, e_init),
            evolv_coeff_β3(mass, n_init, e_init),
            ScaledPeriastronAngle(
                γbar_from_e(e_init),
                γbar2_from_e(e_init, mass),
                γbar3_from_e(e_init, mass),
            ),
        )
    end
end

function n_from_e(coeffs::EvolvCoeffs, e_now::Eccentricity)::MeanMotion
    n0 = coeffs.n_init.n
    e0 = coeffs.e_init.e
    e = e_now.e

    n =
        n0 *
        (e0 / e)^(18 / 19) *
        ((1 - e^2) / (1 - e0^2))^1.5 *
        ((304 + 121 * e0^2) / (304 + 121 * e^2))^(1305 / 2299)
    return MeanMotion(n)
end

function τ_from_t(delay::Time, coeffs::EvolvCoeffs)::ScaledTime
    τ0 = coeffs.τ0.τ
    κ = coeffs.κ
    dt = delay.t
    return ScaledTime(τ0 - κ * dt)
end

function l_from_lbar(
    coeffs::EvolvCoeffs,
    l_init::Angle,
    lbar_now::ScaledMeanAnomaly
)::Angle
    l0 = l_init.θ
    lbar = lbar_now.lbar
    lbar0 = coeffs.lbar0.lbar
    α = coeffs.α
    l = l0 + (lbar0 - lbar) * α
    return Angle(l)
end

function γ_from_γbar(
    coeffs::EvolvCoeffs,
    γ_init::Angle,
    γbar_now::ScaledPeriastronAngle
)::Angle
    γbar1 = γbar_now.γbar1
    γbar2 = γbar_now.γbar2
    γbar3 = γbar_now.γbar3
    γbar10 = coeffs.γbar0_123.γbar1
    γbar20 = coeffs.γbar0_123.γbar2
    γbar30 = coeffs.γbar0_123.γbar3

    β1 = coeffs.β
    β2 = coeffs.β2
    β3 = coeffs.β3

    γ0 = γ_init.θ

    γ = γ0 + (γbar10-γbar1)*β1 + (γbar20-γbar2)*β2 + (γbar30-γbar3)*β3
    return Angle(γ)
end

function evolve_orbit(
    coeffs::EvolvCoeffs,
    l_init::Angle,
    γ_init::Angle,
    delay::Time
)    
    τ::ScaledTime = τ_from_t(delay, coeffs)
    e::Eccentricity = e_from_τ(τ)
    n::MeanMotion = n_from_e(coeffs, e)
    lbar::ScaledMeanAnomaly = lbar_from_e(e)
    γbar123::ScaledPeriastronAngle = γbar123_from_e(e, coeffs.mass)

    l::Angle = l_from_lbar(coeffs, l_init, lbar)
    γ::Angle = γ_from_γbar(coeffs, γ_init, γbar123)

    return n, e, l, γ
end

function derivative_dτ_de(ecc::Eccentricity)::Float64
    e = ecc.e
    return e^(29/19) * (121*e^2+304)^(1181/2299) / (1-e^2)^(3/2)
end

function derivative_de_dt(mass::Mass, norb::MeanMotion, ecc::Eccentricity)::Float64
    Mch = mass.Mch
    n = norb.n
    e = ecc.e
    return (-1/15) * (Mch*n)^(5/3) * n * e * (304+121*e^2) / (1-e^2)^2.5
end

function derivative_dn_dt(mass::Mass, norb::MeanMotion, ecc::Eccentricity)::Float64
    Mch = mass.Mch
    n = norb.n
    e = ecc.e
    return (1/5) * (Mch*n)^(5/3) * n^2 * (96 + 292*e^2 + 37*e^4) / (1-e^2)^3.5
end

function derivative_dlbar_de(ecc::Eccentricity)::Float64
    e = ecc.e
    return e^(11/19) / (304+121*e^2)^(124/2299)
end

function derivative_dγbar_de(ecc::Eccentricity)::Float64
    e = ecc.e
    return e^(-1/19) / (304+121*e^2)^(994/2299)
end

function derivative_dγbar2_de(ecc::Eccentricity, mass::Mass)::Float64
    e = ecc.e
    η = mass.η
    return (e^2*(51-26*η) - 28*η + 78) / e^(13/19) / (304+121*e^2)^(1864/2299)
end

function derivative_dγbar3_de(ecc::Eccentricity, mass::Mass)::Float64
    e = ecc.e
    η = mass.η
    return (18240 − 25376*η + 492*(π^2)*η + 896*η^2 + (28128 − 27840*η + 123*(π^2)*η + 5120*η^2)*e^2 
        + (2496 − 1760*η + 1040*η^2)*e^4 + (1920 − 768*η + (3840 − 1536*η)*e^2)*sqrt(1-e^2)) / e^(25/19) / (304+121*e^2)^(2734/2299)
end

function advance_of_periastron_1PN(mass::Mass, norb::MeanMotion, ecc::Eccentricity)::PeriastronAdvance
    e = ecc.e
    n = norb.n
    M = mass.m
    x = (M*n)^(2/3)
    k = 3*x/(1-e^2)
    return PeriastronAdvance(k)
end

function advance_of_periastron_2PN(mass::Mass, norb::MeanMotion, ecc::Eccentricity)::PeriastronAdvance
    e = ecc.e
    n = norb.n
    M = mass.m
    η = mass.η
    x = (M*n)^(2/3)
    k2 = 0.25 * x^2 * ((51 - 26*η)*e^2 − 28*η + 78) / (1-e^2)^2
    return PeriastronAdvance(k2)
end

function advance_of_periastron_3PN(mass::Mass, norb::MeanMotion, ecc::Eccentricity)::PeriastronAdvance
    e = ecc.e
    n = norb.n
    M = mass.m
    η = mass.η
    x = (M*n)^(2/3)
    k3 = (1/128) * x^3 * (
        (18240 − 25376*η + 492*(π^2)*η + 896*η^2 
        + (28128 − 27840*η + 123*(π^2)*η + 5120*η^2)*e^2 
        + (2496 − 1760*η + 1040*η^2)*e^4 + (1920 − 768*η 
        + (3840 − 1536*η)*e^2)*sqrt(1−e^2))
    ) / (1-e^2)^3
    return PeriastronAdvance(k3)
end