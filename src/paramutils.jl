export mass_from_log10_mass,
    mean_motion_from_log10_freq,
    psrdist_from_pdist,
    Δp_from_deltap,
    dl_from_gwdist,
    mass_from_gwdist

using Unitful
using UnitfulAstro
using PhysicalConstants.CODATA2018

const GMsun = UnitfulAstro.GMsun
const c_0 = CODATA2018.c_0

function mass_from_log10_mass(log10_M::Float64, eta::Float64)::Mass
    M::Float64 = uconvert(u"s", (10.0^log10_M * GMsun) / c_0^3).val
    return Mass(M, eta)
end

function mean_motion_from_log10_freq(log10_F::Float64)::MeanMotion
    return MeanMotion(π * (10.0^log10_F))
end

function psrdist_from_pdist(pdist::Float64)::Distance
    dp = uconvert(u"s", pdist * u"kpc" / c_0).val
    return Distance(dp)
end

function dl_from_gwdist(gwdist::Float64)::Distance
    dp = uconvert(u"s", gwdist * u"Mpc" / c_0).val
    return Distance(dp)
end

function Δp_from_deltap(deltap::Float64)::Time
    year = 365.25 * 24 * 3600
    return Time(-deltap * year)
end

# function mean_motion_from_log10_sidereal_freq(
#     log10_Fb::Float64,
#     e::Eccentricity,
#     mass::Mass,
# )::MeanMotion
#     n0 = MeanMotion(10^log10_Fb)

#     n = n0
#     for _ = 1:3
#         k = advance_of_periastron(mass, n, e)
#         n = MeanMotion(n0.n / (1 + k.k))
#     end

#     return n
# end

function mass_from_gwdist(log10_A, log10_F, e0, gwdist, η)
    n_init = mean_motion_from_log10_freq(log10_F)
    e_init = Eccentricity(e0)

    n0 = n_init.n
    S0 = 10^log10_A
    dl = dl_from_gwdist(gwdist)

    M0 = ((S0 * dl.D / η)^3 * n0)^(1 / 5)
    
    pta_A = m -> gw_amplitude(Mass(m, η), n_init, e_init, dl) / n0
    f = m -> pta_A(m) - S0
    
    ∂S_∂M = (S, M, k) -> (15 + 19*k) * S / (9*(1+k)*M)
    k_aop = m ->  advance_of_periastron(Mass(m, η), n_init, e_init).k
    fp = m -> ∂S_∂M(pta_A(m), m, k_aop(m))

    M = newton_raphson(f, fp, M0)
    return Mass(M, η)
end


function newton_raphson(f, fp, x0, atol=1e-6, maxiter=100)
    for _ in 1:maxiter
        x1 = x0 - f(x0) / fp(x0)

        if abs(x0-x1) < atol
            break
        end

        x0 = x1
    end

    return x0
end