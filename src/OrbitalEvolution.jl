export e_from_tau, tau_from_e

import DataInterpolations
import JLD

function read_precomputed_tau_e(datafile::String)
    data = JLD.load(datafile)
    return data["taus"], data["es"]
end

taus, es = read_precomputed_tau_e("/home/susobhan/Work/GWecc.jl/data/tau_e.jld")
tau_from_e_spline = DataInterpolations.CubicSpline(taus, es)
e_from_tau_spline = DataInterpolations.CubicSpline(es, taus)

emin::Float64, emax::Float64 = minimum(es), maximum(es)
taumin::Float64, taumax::Float64 = minimum(taus), maximum(taus)

a::Float64 = 2 * sqrt(2) / 5 / 5^(63/2299) / 17^(1181/2299);
b::Float64 = 2/sqrt(1-emax) - a*taumax;

function e_from_tau(tau::Float64)
    if tau < 0
        throw(DomainError(tau, "tau < 0 encountered in e_from_tau."));
    elseif tau < taumin
        coeff::Float64 = 2^(559/726) * 3^(19/48) / 19^(145/242);
        e = coeff * tau^(19/48);
    elseif tau > taumax
        e = 1 - 4/(a*tau+b)^2;
    else
        e = e_from_tau_spline(tau);
    end
    
    return e;
end

function tau_from_e(e::Float64)
    if e<0 || e>=1
        throw(DomainError(tau, "e<0 or e>=1 encountered in tau_from_e."));
    elseif e<emin
        coeff = 19 * 19^(1181/2299) / (6*2^(2173/2299));
        tau = coeff * e^(48/19);
    elseif e>emax
        tau = (2/sqrt(1-e) - b)/a;
    else
        tau = tau_from_e_spline(e);
    end

    return tau
end


# fig = plot(taus, es, fmt=:png)
# savefig(fig, "e_vs_tau.png")
