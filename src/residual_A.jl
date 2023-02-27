export residual_PQR, waveform_A, residual_A

"Waveform component functions"
function waveform_A(ecc::Eccentricity, phase::OrbitalPhase)
    e = ecc.e
    su = phase.scu.sinx
    cu = phase.scu.cosx
    s2φ = phase.sc2φ.sinx
    c2φ = phase.sc2φ.cosx

    χ = e * cu
    ξ = e * su

    P = (2 * e^2 - χ^2 + χ - 2) / (1 - χ)^2
    Q = (2 * sqrt(1 - e^2) * ξ) / (1 - χ)^2
    R = χ / (1 - χ)

    hA0 = R
    hA1 = -Q * s2φ + P * c2φ
    hA2 = Q * c2φ + P * s2φ

    return hA0, hA1, hA2
end

"PTA signal component functions."
function residual_PQR(ecc::Eccentricity, scu::SinCos)
    e = ecc.e
    su = scu.sinx
    cu = scu.cosx
    c2u = cu * cu - su * su

    P = (sqrt(1 - e^2) * (c2u - e * cu)) / (1 - e * cu)
    Q = (((e^2 - 2) * cu + e) * su) / (1 - e * cu)
    R = e * su

    return P, Q, R
end

"PTA signal component functions."
function residual_A(ecc::Eccentricity, phase::OrbitalPhase)
    s2ω = phase.sc2ω.sinx
    c2ω = phase.sc2ω.cosx

    P, Q, R = residual_PQR(ecc, phase.scu)

    sA0 = R
    sA1 = -P * s2ω + Q * c2ω
    sA2 = P * c2ω + Q * s2ω

    return sA0, sA1, sA2
end
