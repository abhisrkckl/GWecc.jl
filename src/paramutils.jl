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
