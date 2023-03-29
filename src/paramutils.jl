export mass_from_log10_mass, mean_motion_from_log10_freq, psrdist_from_pdist, Δp_from_deltap

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
