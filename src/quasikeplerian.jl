export advance_of_periastron_1PN,
    advance_of_periastron_2PN, advance_of_periastron_3PN, advance_of_periastron

function advance_of_periastron_1PN(
    mass::Mass,
    norb::MeanMotion,
    ecc::Eccentricity,
)::PeriastronAdvance
    e = ecc.e
    n = norb.n
    M = mass.m
    x = (M * n)^(2 / 3)
    k = 3 * x / (1 - e^2)
    return PeriastronAdvance(k)
end

function advance_of_periastron_2PN(
    mass::Mass,
    norb::MeanMotion,
    ecc::Eccentricity,
)::PeriastronAdvance
    e = ecc.e
    n = norb.n
    M = mass.m
    η = mass.η
    x = (M * n)^(2 / 3)
    k2 = 0.25 * x^2 * ((51 - 26 * η) * e^2 − 28 * η + 78) / (1 - e^2)^2
    return PeriastronAdvance(k2)
end

function advance_of_periastron_3PN(
    mass::Mass,
    norb::MeanMotion,
    ecc::Eccentricity,
)::PeriastronAdvance
    e = ecc.e
    n = norb.n
    M = mass.m
    η = mass.η
    x = (M * n)^(2 / 3)
    k3 =
        (1 / 128) *
        x^3 *
        ((
            18240 − 25376 * η +
            492 * (π^2) * η +
            896 * η^2 +
            (28128 − 27840 * η + 123 * (π^2) * η + 5120 * η^2) * e^2 +
            (2496 − 1760 * η + 1040 * η^2) * e^4 +
            (1920 − 768 * η + (3840 − 1536 * η) * e^2) * sqrt(1 − e^2)
        )) / (1 - e^2)^3
    return PeriastronAdvance(k3)
end

function advance_of_periastron(mass::Mass, n::MeanMotion, ecc::Eccentricity)
    k1 = advance_of_periastron_1PN(mass, n, ecc).k
    k2 = advance_of_periastron_2PN(mass, n, ecc).k
    k3 = advance_of_periastron_3PN(mass, n, ecc).k
    return PeriastronAdvance(k1 + k2 + k3)
end
