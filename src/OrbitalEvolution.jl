export e_from_tau, tau_from_e

import DataInterpolations
import JLD
using HypergeometricFunctions

function read_precomputed_tau_e(datafile::String)
    data = JLD.load(datafile)
    return data["taus"], data["es"]
end

taus, es = read_precomputed_tau_e("/home/susobhan/Work/GWecc.jl/data/tau_e.jld")
tau_from_e_spline = DataInterpolations.CubicSpline(taus, es)
e_from_tau_spline = DataInterpolations.CubicSpline(es, taus)

struct ScaledTime
    τ::Float64
    ScaledTime(τ::Float64) = τ>=0 ? new(τ) : throw(DomainError(τ, "τ<0 encountered."))
end

struct Eccentricity
    e::Float64
    Eccentricity(e::Float64) = (e>0 && e<1) ? new(e) : throw(DomainError(e, "e out of range."))
end

struct Mass
    m::Float64
    Mass(m::Float64) = (m>=5e-2 && m<=5e4) ? new(m) : throw(DomainError(m, "m out of range."))
end

struct MeanMotion
    n::Float64
    MeanMotion(n::Float64) = (n>=6.3e-10 && n<=6.3e-6) ? new(n) : throw(DomainError(n, "n out of range."))
end

struct ScaledAngle
    θbar::Float64
    ScaledAngle(θbar::Float64) = (θbar>=0 && θbar<=0.5) ? new(θbar) : throw(DomainError(θbar, "θbar out of range."))
end

struct SymmetricMassRatio
    η::Float64
    ScaledAngle(η::Float64) = (η>0 && η<=0.25) ? new(η) : throw(DomainError(η, "η out of range."))
end

eccmin = Eccentricity(minimum(es)) 
eccmax = Eccentricity(maximum(es))
taumin = ScaledTime(minimum(taus))
taumax = ScaledTime(maximum(taus))

a::Float64 = 2 * sqrt(2) / 5 / 5^(63/2299) / 17^(1181/2299);
b::Float64 = 2/sqrt(1-eccmax.e) - a*taumax.τ;

macro validate(value, condition, message)
    if not condition(value)
        throw(DomainError(value, message))
    end
end

function e_from_τ(tau::ScaledTime) :: Eccentricity
    τ = tau.τ
    τmin = taumin.τ
    τmax = taumax.τ

    if τ < τmin
        coeff::Float64 = 2^(559/726) * 3^(19/48) / 19^(145/242);
        e = coeff * τmin^(19/48);
    elseif τ > τmax
        e = 1 - 4/(a*τ+b)^2;
    else
        e = e_from_tau_spline(τ);
    end
    
    return Eccentricity(e);
end

function tau_from_e(ecc::Eccentricity) :: ScaledTime
    e = ecc.e
    emin = eccmin.e
    emax = eccmax.e
    
    if e<emin
        coeff = 19 * 19^(1181/2299) / (6*2^(2173/2299));
        τ = coeff * e^(48/19);
    elseif e>emax
        τ = (2/sqrt(1-e) - b)/a;
    else
        τ = tau_from_e_spline(e);
    end

    return ScaledTime(τ)
end

function evolv_coeff_κ(Mchirp::Mass, n_init::MeanMotion, e_init::Eccentricity) :: Float64
    Mch = Mchirp.m
    n0 = n_init.n
    e0 = e_init.e
    
    return Mch^(5/3) / 15 * n0^(8/3) * e0^(48/19) * (304 + 121*e0^2)^(3480/2299) / (1-e0^2)^4;   
end

function evolv_coeff_α(Mchirp::Mass, n_init::MeanMotion, e_init::Eccentricity) :: Float64
    Mch = Mchirp.m
    n0 = n_init.n
    e0 = e_init.e

    return 15 / Mchirp^(5/3) / n0^(5/3) / e0^(30/19) / (304 + 121*e0^2)^(2175/2299) * (1-e0^2)^2.5;
end

function evolv_coeff_β(Mchirp::Mass, Mtot::Mass, n_init::MeanMotion, e_init::Eccentricity) :: Float64
    Mch = Mchirp.m
    M = Mtot.m
    n0 = n_init.n
    e0 = e_init.e

    return 45 / Mch^(5/3) * M^(2/3) / n0 / e0^(18/19) / pow(304+121*e0^2)^(1305/2299) * (1-e0^2)^1.5; 
end

function evolv_coeff_β2(Mchirp::Mass, Mtot::Mass, n_init::MeanMotion, e_init::Eccentricity) :: Float64
    Mch = Mchirp.m
    M = Mtot.m
    n0 = n_init.n
    e0 = e_init.e
    
    return 15 / (4*Mch^(5/3)) * M^(4/3) / cbrt(n0) / e0^(6/19) / (304+121*e0^2)^(435/2299) * sqrt(1-e0^2);
end

function evolv_coeff_β3(Mchirp::Mass, Mtot::Mass, n_init::MeanMotion, e_init::Eccentricity) :: Float64
    Mch = Mchirp.m
    M = Mtot.m
    n0 = n_init.n
    e0 = e_init.e
    
    return 15 / (128*Mchirp^(5/3)) * (M*M) * cbrt(n0) * e0^(6/19) * (304+121*e0^2)^(435/2299) / sqrt(1-e0^2);
end

function n_from_e(n_init::MeanMotion, e_init::Eccentricity, e_now::Eccentricity) :: MeanMotion
    n0 = n_init.n
    e0 = e_init.e
    e = e_now.e
    
    n = n0 * (e0/e)^(18/19) * ((1-e2)/(1-e02))^1.5 * ((304 + 121*e02)/(304 + 121*e2))^(1305/2299);
    return MeanMotion(n)
end

function lbar_from_e(ecc::Eccentricity) :: ScaledAngle
    e = ecc.e

    coeff::Float64 = 19^(2175/2299) / 30 / 2^(496/2299);
    lbar = coeff * e * e^(11/19) * _₂F₁(124/2299, 15/19, 34/19, -121*e^2/304);
    return ScaledAngle(lbar)
end

function γbar_from_e(ecc::Eccentricity) :: ScaledAngle
    e = ecc.e
    coeff::Float64 = 19^(1305/2299) / 36 / 2^(1677/2299);
    γbar = coeff * e^(18/19) * _₂F₁(994/2299,  9/19, 28/19, -121*e^2/304);
    return ScaledAngle(γbar)
end

function γbar2_from_e(ecc::Eccentricity, eta::SymmetricMassRatio) :: ScaledAngle
    e = ecc.e
    η = eta.η
    coeff::Float64 = 3 * 2^(1740/2299) * 19^(435/2299);
    γbar2 = e^(6/19) / 336 * (  4 * (51-26*η) * (304+121*e^2)^(435/2299)
                                + coeff * (23+2*eta) * _₂F₁(3/19, 1864/2299, 22/19, -121*e^2/304)
                            )
end

function γbar3_from_e(ecc::Eccentricity, eta::SymmetricMassRatio) :: ScaledAngle
    """
    γbar3(e) = (γbar3_0 + η*γbar3_1(e) + (η^2)*γbar3_2(e)) / e^(6/19)
     
    The functions gbar30 and gbar31 contain Appell F1 functions that are
    hard to compute. So I am using Pade approximants here. Since this is a
    3PN correction the error due to this approximation should be negligible.
    This probably won't work if e == 0.
    """
    e = ecc.e
    η = eta.η
    γbar3_0::Float64 = -71.19148913450734 + (14.212810668086124*(e*e) - 2.7809293422212327*(e^4) - 0.5247248712090506*(e^6))/(1. + 0.021216346439953203*(e*e) - 0.08582538305359699*(e^4));
    γbar3_1::Float64 = 75.17536852020817 + (-11.720605504532221*(e*e) + 0.9937219102426796*(e^4) + 0.3416818503954963*(e^6))/(1. + 0.11825867242669967*(e*e) - 0.05881762636038525*(e^4));
    γbar3_2::Float64 = -3.1640661837558817 + (3.109257298381722*(e*e) + 1.0709804529990066*(e^4) + 0.05588268809419086*(e^6))/(1. + 0.46114315553206664*(e*e) + 0.041157522000398405*(e^4));
    γbar3 = (γbar3_0 + η*γbar3_1 * η*η*γbar3_2) / e^(6/19);
    return ScaledAngle(γbar3)
end

