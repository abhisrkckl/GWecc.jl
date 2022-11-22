export e_from_τ, τ_from_e, eccmin, eccmax, taumin, taumax
export evolv_coeff_κ,
    evolv_coeff_α, evolv_coeff_β, evolv_coeff_β2, evolv_coeff_β3, EvolvCoeffs
export n_from_e, lbar_from_e, γbar_from_e, γbar2_from_e, γbar3_from_e

import DataInterpolations
import JLD
using HypergeometricFunctions

include("Parameters.jl")

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

function evolv_coeff_κ(mass::Mass, n_init::MeanMotion, e_init::Eccentricity)::Float64
    Mch = mass.Mch
    n0 = n_init.n
    e0 = e_init.e

    return Mch^(5 / 3) / 15 * n0^(8 / 3) * e0^(48 / 19) * (304 + 121 * e0^2)^(3480 / 2299) /
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

function n_from_e(n_init::MeanMotion, e_init::Eccentricity, e_now::Eccentricity)::MeanMotion
    n0 = n_init.n
    e0 = e_init.e
    e = e_now.e

    n =
        n0 *
        (e0 / e)^(18 / 19) *
        ((1 - e^2) / (1 - e0^2))^1.5 *
        ((304 + 121 * e0^2) / (304 + 121 * e^2))^(1305 / 2299)
    return MeanMotion(n)
end

function lbar_from_e(ecc::Eccentricity)::ScaledMeanAnomaly
    e = ecc.e

    coeff::Float64 = 19^(2175 / 2299) / 30 / 2^(496 / 2299)
    lbar = coeff * e * e^(11 / 19) * _₂F₁(124 / 2299, 15 / 19, 34 / 19, -121 * e^2 / 304)
    return ScaledMeanAnomaly(lbar)
end

function γbar_from_e(ecc::Eccentricity)::Float64
    e = ecc.e
    coeff::Float64 = 19^(1305 / 2299) / 36 / 2^(1677 / 2299)
    γbar = coeff * e^(18 / 19) * _₂F₁(994 / 2299, 9 / 19, 28 / 19, -121 * e^2 / 304)
    return γbar
end

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

function τ_from_t(delay::Time, tau0::ScaledTime, κ::Float64)::ScaledTime
    return ScaledTime(tau0.τ - κ * delay.t)
end

function l_from_lbar(
    lbar_now::ScaledMeanAnomaly,
    lbar_init::ScaledMeanAnomaly,
    coeff::Float64,
    l_init::Angle,
)::Angle
    l0 = l_init.θ
    lbar = lbar_now.lbar
    lbar0 = lbar_init.lbar
    l = l0 + (lbar0 - lbar) * coeff
    return Angle(l)
end

function γ_from_γbar(
    γbar_now::ScaledPeriastronAngle,
    γbar_init::ScaledPeriastronAngle,
    coeffs::Float64,
)::Angle
    γbar1 = γbar_now.γbar1
    γbar2 = γbar_now.γbar2
    γbar3 = γbar_now.γbar3
    γbar10 = γbar_init.γbar1
    γbar20 = γbar_init.γbar2
    γbar30 = γbar_init.γbar3

    θbar = θbar_now.θbar
    θbar0 = θbar_init.θbar
    θ = (θbar0 - θbar) * coeff
    return Angle(θ)
end

struct EvolvCoeffs
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

#= function evolve_orbit(
    Mtot::Mass,
    eta::SymmetricMassRatio,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l_init::Angle,
    proj_pars::ProjectionParams,
    delay::Time,
    pulsar_term::Bool,
)

    γ_init = Angle(pulsar_term ? proj_pars.γp : proj_pars.γ0)

    Mchirp::Mass = chirp_mass(Mtot, eta)

    κ::Float64 = evolv_coeff_κ(Mchirp, n_init, e_init)
    tau0::ScaledTime = τ_from_e(e_init)
    tau::ScaledTime = τ_from_t(delay, tau0, κ)

    if tau.τ < 0
        throw(DomainError(tau, "The binary has already merged (τ<0)."))
    end

    e::Eccentricity = e_from_τ(tau)
    n::MeanMotion = n_from_e(n_init, e_init, e)

    α::Float64 = evolv_coeff_α(Mchirp, n_init, e_init)
    lbar0::ScaledAngle = lbar_from_e(e_init)
    lbar::ScaledAngle = lbar_from_e(e)
    l::Angle = θ_from_θbar(lbar, lbar0, α, l_init)

    β::Float64 = evolv_coeff_β(Mchirp, Mtot, n_init, e_init)
    γbar0::ScaledAngle = γbar_from_e(e0)
    γbar::ScaledAngle = γbar_from_e(e)
    γ1::Angle = θ_from_θbar(γbar, γbar0, β, γ_init)

    β2::Float64 = evolv_coeff_β2(Mchirp, Mtot, n_init, e_init)
    γbar20::ScaledAngle = γbar2_from_e(e_init, eta)
    γbar2::ScaledAngle = γbar2_from_e(e, eta)
    γ2::Angle = θ_from_θbar(γbar2, γbar20, β2)

    β3::Float64 = evolv_coeff_β3(Mchirp, Mtot, n_init, e_init)
    γbar30::ScaledAngle = γbar3_from_e(e_init, eta)
    γbar3::ScaledAngle = γbar3_from_e(e, eta)
    γ3::Angle = θ_from_θbar(γbar3, γbar30, β3)

    γ = Angle(γ1.θ + γ2.θ + γ3.θ)

    return n, e, l, γ
end
 =#
