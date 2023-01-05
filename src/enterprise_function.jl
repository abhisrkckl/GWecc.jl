function eccentric_pta_signal_planck18_1psr(
    toas,
    pdist,
    alpha,
    psi,
    cos_inc,
    log10_M,
    eta,
    log10_F,
    e0,
    gamma0,
    gammap,
    l0,
    lp,
    tref,
    log10_zc,
    psrTerm::Bool=false
)
    mass = mass_from_log10_mass(log10_M, eta)
    n_init = mean_motion_from_log10_freq(π * (10.0^log10_F))
    e_init = Eccentricity(e0)
    l0p = InitPhaseParams(l0, lp)
    proj = ProjectionParams(psi, cos_inc, gamma0, gammap)
    z, dl = redshift_luminosity_dist_from_log10_redshift(log10_zc)
    dp = psrdist_from_pdist(pdist)
    α = AzimuthParam(alpha)
    terms = psrTerm ? [EARTH, PULSAR] : [EARTH]
    tref = Time(tref)
    tEs = Time.(toas)

    return residuals_1psr(
        mass,
        n_init,
        e_init,
        l0p,
        proj,
        dl,
        dp,
        α,
        z,
        terms,
        tref,
        tEs
    )
end