export eccentric_pta_signal_planck18_1psr

function eccentric_pta_signal_planck18_1psr(
    toas,
    pdist::Float64,
    alpha::Float64,
    psi::Float64,
    cos_inc::Float64,
    log10_M::Float64,
    eta::Float64,
    log10_F::Float64,
    e0::Float64,
    gamma0::Float64,
    gammap::Float64,
    l0::Float64,
    lp::Float64,
    tref::Float64,
    log10_zc::Float64,
    psrTerm::Bool = false,
)
    mass = mass_from_log10_mass(log10_M, eta)
    n_init = mean_motion_from_log10_freq(log10_F)
    e_init = Eccentricity(e0)
    l0p = InitPhaseParams(l0, lp)
    proj = ProjectionParams(psi, cos_inc, gamma0, gammap)
    z, dl = redshift_luminosity_dist_from_log10_redshift(log10_zc)
    dp = psrdist_from_pdist(pdist)
    α = AzimuthParam(alpha)
    terms = psrTerm ? [EARTH, PULSAR] : [EARTH]
    tref = Time(tref)
    tEs = Time.(toas)

    return residuals_1psr(mass, n_init, e_init, l0p, proj, dl, dp, α, z, terms, tref, tEs)
end
