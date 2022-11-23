export derivative_dτ_de, derivative_de_dt, derivative_dn_dt
export derivative_dlbar_de, derivative_dγbar_de, derivative_dγbar2_de, derivative_dγbar3_de

function derivative_dτ_de(ecc::Eccentricity)::Float64
    e = ecc.e
    return e^(29 / 19) * (121 * e^2 + 304)^(1181 / 2299) / (1 - e^2)^(3 / 2)
end

function derivative_de_dt(mass::Mass, norb::MeanMotion, ecc::Eccentricity)::Float64
    Mch = mass.Mch
    n = norb.n
    e = ecc.e
    return (-1 / 15) * (Mch * n)^(5 / 3) * n * e * (304 + 121 * e^2) / (1 - e^2)^2.5
end

function derivative_dn_dt(mass::Mass, norb::MeanMotion, ecc::Eccentricity)::Float64
    Mch = mass.Mch
    n = norb.n
    e = ecc.e
    return (1 / 5) * (Mch * n)^(5 / 3) * n^2 * (96 + 292 * e^2 + 37 * e^4) / (1 - e^2)^3.5
end

function derivative_dlbar_de(ecc::Eccentricity)::Float64
    e = ecc.e
    return e^(11 / 19) / (304 + 121 * e^2)^(124 / 2299)
end

function derivative_dγbar_de(ecc::Eccentricity)::Float64
    e = ecc.e
    return e^(-1 / 19) / (304 + 121 * e^2)^(994 / 2299)
end

function derivative_dγbar2_de(ecc::Eccentricity, mass::Mass)::Float64
    e = ecc.e
    η = mass.η
    return (e^2 * (51 - 26 * η) - 28 * η + 78) / e^(13 / 19) /
           (304 + 121 * e^2)^(1864 / 2299)
end

function derivative_dγbar3_de(ecc::Eccentricity, mass::Mass)::Float64
    e = ecc.e
    η = mass.η
    return (
        18240 − 25376 * η +
        492 * (π^2) * η +
        896 * η^2 +
        (28128 − 27840 * η + 123 * (π^2) * η + 5120 * η^2) * e^2 +
        (2496 − 1760 * η + 1040 * η^2) * e^4 +
        (1920 − 768 * η + (3840 − 1536 * η) * e^2) * sqrt(1 - e^2)
    ) / e^(25 / 19) / (304 + 121 * e^2)^(2734 / 2299)
end
