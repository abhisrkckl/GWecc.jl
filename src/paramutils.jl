using Cosmology
using Unitful
using UnitfulAstro
using PhysicalConstants.CODATA2018

GMsun = UnitfulAstro.GMsun
c_0 = CODATA2018.c_0

# Based on astropy.cosmology.Planck18
planck18 = cosmology(
    h = 0.6766,
    Neff = 3.046,
    OmegaK = 0.0,
    OmegaM = 0.30966,
    OmegaR = 5.402015137139352e-5,
    Tcmb = 2.7255,
    w0 = -1,
    wa = 0
)

function mass_from_log10_mass(log10_M::Float64, eta::Float64)::Mass
    M = uconvert(u"s", (10.0^log10_M * GMsun) / c_0^3).val
    return Mass(M, eta)
end

function mean_motion_from_log10_freq(log10_F::Float64)::MeanMotion
    return MeanMotion(π * (10.0^log10_F))
end

function redshift_luminosity_dist_from_log10_redshift(log10_z::Float64)::Tuple{Redshift,Distance}
    z = 10^log10_z
    dl = uconvert(u"s", luminosity_dist(planck18, z) / c_0).val
    return Redshift(z), Distance(dl)
end

function psrdist_from_pdist(pdist::Float64)::Distance
    dp = uconvert(u"s", pdist * u"kpc" / c_0).val
    return Distance(dp)
end