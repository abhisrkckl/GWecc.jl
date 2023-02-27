export gw_amplitude, gwres_amplitude_ratio, waveform_coeffs_c

"GW amplitude"
function gw_amplitude(
    mass::Mass,
    norb::MeanMotion,
    ecc::Eccentricity,
    dl::Distance,
)::Float64
    m, η = mass.m, mass.η
    dgw = dl.D
    x = pn_param_x(mass, norb, ecc).x
    return m * η * x / dgw
end

# function gw_amplitude_ratio(
#     mass::Mass,
#     n0::MeanMotion,
#     e0::Eccentricity,
#     n1::MeanMotion,
#     e1::Eccentricity,
# )
#     x0 = pn_param_x(mass, n0, e0).x
#     x1 = pn_param_x(mass, n1, e1).x
#     return x1 / x0
# end

function gwres_amplitude_ratio(
    mass::Mass,
    n0::MeanMotion,
    e0::Eccentricity,
    n1::MeanMotion,
    e1::Eccentricity,
)
    x0 = pn_param_x(mass, n0, e0).x
    x1 = pn_param_x(mass, n1, e1).x
    return (x1 / x0) / (n1.n / n0.n)
end

"Waveform coefficients that depend on the inclination"
function waveform_coeffs_c(proj::ProjectionParams)
    ci = proj.cosι
    return 1 - ci^2, 1 + ci^2, 2 * ci
end
