export gw_amplitude, residual

function gw_amplitude(mass::Mass, norb::MeanMotion, ecc::Eccentricity, dl::Distance)
    m, η = mass.m, mass.η
    dgw = dl.D
    x = pn_param_x(mass, norb, ecc).x
    return m*η*x/dgw;
end

function residual_A(ecc::Eccentricity, phase::OrbitalPhase)
    e = ecc.e
    su = phase.scu.sinx
    cu = phase.scu.cosx
    c2u = cu*cu - su*su
    s2ω = phase.sc2ω.sinx
    c2ω = phase.sc2ω.cosx

    P = ((e + (-2 + e*e)*cu)*su)/(1 - e*cu)
    Q = (sqrt(1-e^2)*(e*cu - c2u))/(1 - e*cu)
    R = e*su

    A0 = R
    A1 = -P*s2ω + Q*c2ω
    A2 =  P*c2ω + Q*s2ω

    return A0, A1, A2
end

function residual_a_coeffs(proj::ProjectionParams)
    ci = proj.cosι
    return 1-ci^2, 1+ci^2, 2*ci
end

function residual_spx(mass::Mass, coeffs::EvolvCoeffs, l_init::Angle, proj::ProjectionParams, dl::Distance, dt::Time, psrterm::Bool)
    γ_init = psrterm ? Angle(proj.γp) : Angle(proj.γ0)
    
    n, e, l, γ = evolve_orbit(coeffs, l_init, γ_init, dt)
    phase = OrbitalPhase(mass, n, e, l, γ)

    A0, A1, A2 = residual_A(e, phase)
    a0, a1, a2 = residual_a_coeffs(proj)

    h0 = gw_amplitude(mass, n, e, dl)
    s0 = h0 / n.n
    
    sA = s0 * (a1*A1 + a0*A0)
    sB = s0 * (a2*A2)
    
    s2ψ, c2ψ = proj.sc2ψ.sinx, proj.sc2ψ.cosx
    
    sp = c2ψ*sA - s2ψ*sB
    sx = s2ψ*sA + c2ψ*sB

    return sp, sx
end

function residual(mass::Mass, coeffs::EvolvCoeffs, l_init::Angle, proj::ProjectionParams, dl::Distance, dt::Time, ap::AntennaPattern, terms::Vector{Term}, Δp::Time)
    sp = 0.0
    sx = 0.0
    
    if EARTH in terms
        spE, sxE = residual_spx(mass, coeffs, l_init, proj, dl, dt, false)
        sp = sp - spE
        sx = sx - sxE
    end

    if PULSAR in terms
        dtp = Time(dt.t + Δp.t)
        spP, sxP = residual_spx(mass, coeffs, l_init, proj, dl, dtp, true)
        sp = sp + spP
        sx = sx + sxP
    end

    return ap.Fp*sp + ap.Fx*sx
end
