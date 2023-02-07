export waveform_PQR, residual_PQR

"Waveform component functions (PQR)."
function waveform_PQR(ecc::Eccentricity, scu::SinCos)
    e = ecc.e
    su = scu.sinx
    cu = scu.cosx
    
    sf = sqrt(1-e^2)*su/(1-e*cu)
    cf = (cu-e)/(1-e*cu)
    s2f = 2*sf*cf
    c2f = cf*cf - sf*sf

    χ = e * cu
    ξ = e * su

    p = (2 * e^2 - χ^2 + χ - 2) / (1 - χ)^2
    q = (2 * sqrt(1 - e^2) * ξ) / (1 - χ)^2
    r = χ / (1 - χ)

    P = p*s2f + q*c2f
    Q = p*c2f - q*s2f
    R = r

    return P, Q, R
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