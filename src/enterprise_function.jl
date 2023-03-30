export eccentric_pta_signal, eccentric_pta_signal_1psr, eccentric_pta_signal_target, validate_params, validate_params_target

"High-level interface for the residuals function. This function computes 
the eccentric PTA signal for a single pulsar given a collection of TOAs 
and source parameters.

Parameters:
    toas : Collection of TOAs (s)
    theta : Zenith angle of the pulsar (rad)
    phi : RA of the pulsar (rad)
    pdist : Pulsar distance (kpc)
    cos_gwtheta : Cos zenith angle of the GW source
    gwphi : RA of the GW source (rad)
    psi : GW polarization angle (rad)
    cos_inc : Cos inclination
    log10_M : Log10 total mass (Msun)
    eta : Symmetric mass ratio
    log10_F : Log10 initial GW frequency (Hz)
    e0 : Initial eccentricity
    gamma0 : Initial periastron angle for the Earth term (rad)
    gammap : Initial periastron angle for the Pulsar term (rad)
    l0 : Initial mean anomaly for the Earth term (rad)
    lp : Initial mean anomaly for the Pulsar term (rad)
    tref : Fiducial time (s)
    log10_A : Log10 PTA signal amplitude (s)
    psrTerm : Whether to include the pulsar term
    spline : Whether to use spline-based fast computation
"
function eccentric_pta_signal(
    toas,
    theta::Float64,
    phi::Float64,
    pdist::Float64,
    cos_gwtheta::Float64,
    gwphi::Float64,
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
    log10_A::Float64,
    psrTerm::Bool = false,
    spline::Bool = false,
)
    mass = mass_from_log10_mass(log10_M, eta)

    n_init = mean_motion_from_log10_freq(log10_F)
    e_init = Eccentricity(e0)
    # n_init = mean_motion_from_log10_sidereal_freq(log10_F, e_init, mass)

    l0p = InitPhaseParams(l0, lp)

    proj = ProjectionParams(10^log10_A, psi, cos_inc, gamma0, gammap)

    dp = psrdist_from_pdist(pdist)

    ra_psr = phi
    dec_psr = π / 2 - theta
    psrpos = SkyLocation(ra_psr, dec_psr)

    ra_gw = gwphi
    dec_gw = asin(cos_gwtheta)
    gwpos = SkyLocation(ra_gw, dec_gw)

    terms = psrTerm ? [EARTH, PULSAR] : [EARTH]
    tref = Time(tref)
    tEs = Time.(toas)

    res = spline ? residuals_spline : residuals

    return res(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, terms, tref, tEs)
end

function eccentric_pta_signal_target(
    toas,
    theta::Float64,
    phi::Float64,
    pdist::Float64,
    cos_gwtheta::Float64,
    gwphi::Float64,
    psi::Float64,
    cos_inc::Float64,
    eta::Float64,
    log10_F::Float64,
    e0::Float64,
    gamma0::Float64,
    gammap::Float64,
    l0::Float64,
    lp::Float64,
    tref::Float64,
    log10_A::Float64,
    gwdist::Float64,
    psrTerm::Bool = false,
    spline::Bool = false,
)
    n_init = mean_motion_from_log10_freq(log10_F)
    e_init = Eccentricity(e0)

    mass = mass_from_gwdist(log10_A, log10_F, e0, gwdist, eta)

    # n_init = mean_motion_from_log10_sidereal_freq(log10_F, e_init, mass)

    l0p = InitPhaseParams(l0, lp)

    proj = ProjectionParams(10^log10_A, psi, cos_inc, gamma0, gammap)

    dp = psrdist_from_pdist(pdist)

    ra_psr = phi
    dec_psr = π / 2 - theta
    psrpos = SkyLocation(ra_psr, dec_psr)

    ra_gw = gwphi
    dec_gw = asin(cos_gwtheta)
    gwpos = SkyLocation(ra_gw, dec_gw)

    terms = psrTerm ? [EARTH, PULSAR] : [EARTH]
    tref = Time(tref)
    tEs = Time.(toas)

    res = spline ? residuals_spline : residuals

    return res(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, terms, tref, tEs)
end

"High-level interface for the residuals_1psr function. This function computes 
the eccentric PTA signal for a single pulsar given a collection of TOAs and 
source parameters.

Parameters:
    toas : Collection of TOAs (s)
    sigma : Projection angle 1 (rad)
    rho : Projection angle 2 (rad)
    log10_M : Log10 total mass (Msun)
    eta : Symmetric mass ratio
    log10_F : Log10 initial GW frequency (Hz)
    e0 : Initial Eccentricity
    l0 : Initial mean anomaly (rad)
    tref : Fiducial time (s)
    log10_A : Log10 PTA signal amplitude (s)
    deltap : Pulsar term delay (yr)
    psrTerm : Whether to include the pulsar term
    spline : Whether to use spline-based fast computation
"
function eccentric_pta_signal_1psr(
    toas,
    sigma::Float64,
    rho::Float64,
    log10_M::Float64,
    eta::Float64,
    log10_F::Float64,
    e0::Float64,
    l0::Float64,
    tref::Float64,
    log10_A::Float64,
    deltap::Float64,
    psrTerm::Bool = false,
    spline::Bool = false,
)
    mass = mass_from_log10_mass(log10_M, eta)
    n_init = mean_motion_from_log10_freq(log10_F)
    e_init = Eccentricity(e0)
    # n_init = mean_motion_from_log10_sidereal_freq(log10_F, e_init, mass)
    l_init = Angle(l0)
    proj = ProjectionParams1psr(10^log10_A, sigma, rho)
    Δp = Δp_from_deltap(deltap)
    terms = psrTerm ? [EARTH, PULSAR] : [EARTH]
    tref = Time(tref)
    tEs = Time.(toas)

    res = spline ? residuals_1psr_spline : residuals_1psr

    return res(mass, n_init, e_init, l_init, proj, Δp, terms, tref, tEs)
end

function validate_params(
    log10_M::Float64,
    eta::Float64,
    log10_F::Float64,
    e0::Float64,
    tref::Float64,
    tmax::Float64,
)
    try
        delay = Time(tmax - tref)
        mass = mass_from_log10_mass(log10_M, eta)
        n_init = mean_motion_from_log10_freq(log10_F)
        e_init = Eccentricity(e0)
        l_init = Angle(0.0)
        γ_init = Angle(0.0)
        coeffs = EvolvCoeffs(mass, n_init, e_init)

        # Throws if binary merges
        n, e, l, γ = evolve_orbit(coeffs, l_init, γ_init, delay)

        # Throws if QKP fails
        phase = OrbitalPhase(mass, n, e, l, γ)
    catch ex
        return false
    end

    return true
end

function validate_params_target(
    log10_A::Float64,
    eta::Float64,
    log10_F::Float64,
    e0::Float64,
    gwdist::Float64,
    tref::Float64,
    tmax::Float64,
)
    try
        delay = Time(tmax - tref)
        mass = mass_from_gwdist(log10_A, log10_F, e0, gwdist, eta)
        n_init = mean_motion_from_log10_freq(log10_F)
        e_init = Eccentricity(e0)
        l_init = Angle(0.0)
        γ_init = Angle(0.0)
        coeffs = EvolvCoeffs(mass, n_init, e_init)

        # Throws if binary merges
        n, e, l, γ = evolve_orbit(coeffs, l_init, γ_init, delay)

        # Throws if QKP fails
        phase = OrbitalPhase(mass, n, e, l, γ)
    catch ex
        return false
    end

    return true
end
