export PNParam,
    pn_param_x,
    PeriastronAdvance,
    advance_of_periastron,
    advance_of_periastron_1PN,
    advance_of_periastron_2PN,
    advance_of_periastron_3PN,
    angular_eccentricity

"PN expansion parameter"
struct PNParam
    x::Float64
    PNParam(x::Float64) =
        (x > 0 && x < 0.25) ? new(x) : throw(DomainError(x, "x out of range."))
end

"Periastron advance per orbit"
struct PeriastronAdvance
    k::Float64
    PeriastronAdvance(k::Float64) =
        abs(k) < 0.25 ? new(k) : throw(DomainError(k, "k out of range."))
end

"PN expansion parameter"
function pn_param_x(mass::Mass, norb::MeanMotion, ecc::Eccentricity)::PNParam
    m = mass.m
    n = norb.n
    e = ecc.e
    x = (m * n)^(2 / 3)
    k = advance_of_periastron(mass, norb, ecc).k
    x = x * (1 + k)^(2 / 3)
    return PNParam(x)
end

"Advance of periastron at 1PN order"
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

"Advance of periastron at 2PN order"
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

"Advance of periastron at 3PN order"
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

"Advance of periastron accurate up to 3PN order"
function advance_of_periastron(mass::Mass, n::MeanMotion, ecc::Eccentricity)
    k1 = advance_of_periastron_1PN(mass, n, ecc).k
    k2 = advance_of_periastron_2PN(mass, n, ecc).k
    k3 = advance_of_periastron_3PN(mass, n, ecc).k
    return PeriastronAdvance(k1 + k2 + k3)
end

"Angular eccentricity accurate up to 2PN order"
function angular_eccentricity(mass::Mass, n::MeanMotion, ecc::Eccentricity)
    x = pn_param_x(mass, n, ecc).x
    e = ecc.e
    η = mass.η
    ots = sqrt(1 - e^2)
    eφ =
        e * (
            1 +
            x * (4 - η) +
            (x^2) * (
                4 * (-12 * (26 + 15 * ots) + η * (17 + 72 * ots + η)) +
                (e^2) * (1152 + η * (-659 + 41 * η))
            ) / (96 * (-1 + e^2))
        )
    return Eccentricity(eφ)
end
