export spline_time_samples, residuals_spline, residuals_1psr_spline

using CubicHermiteSpline

function spline_time_samples(tEs::Vector{Time})
    day = 24.0 * 3600.0
    tEs_day = [t.t / day for t in tEs]
    append!(tEs_day, [minimum(tEs_day) - 1, maximum(tEs_day) + 1])
    return Time.(sort(unique(round.(tEs_day))) * day)
end

function residuals_spline(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    dl::Distance,
    dp::Distance,
    psrpos::SkyLocation,
    gwpos::SkyLocation,
    z::Redshift,
    terms::Vector{Term},
    tref::Time,
    tspl::Vector{Time},
    tEs::Vector{Time},
)
    ss, hs = residuals_and_waveform(
        mass,
        n_init,
        e_init,
        l0p,
        proj,
        dl,
        dp,
        psrpos,
        gwpos,
        z,
        terms,
        tref,
        tspl,
    )
    spl = CubicHermiteSplineInterpolation([t.t for t in tspl], ss, hs)
    return interp(spl, [t.t for t in tEs])
end

function residuals_spline(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    dl::Distance,
    dp::Distance,
    psrpos::SkyLocation,
    gwpos::SkyLocation,
    z::Redshift,
    terms::Vector{Term},
    tref::Time,
    tEs::Vector{Time},
)
    tspl = spline_time_samples(tEs)
    return residuals_spline(
        mass,
        n_init,
        e_init,
        l0p,
        proj,
        dl,
        dp,
        psrpos,
        gwpos,
        z,
        terms,
        tref,
        tspl,
        tEs,
    )
end

function residuals_1psr_spline(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    dl::Distance,
    dp::Distance,
    α::AzimuthParam,
    z::Redshift,
    terms::Vector{Term},
    tref::Time,
    tspl::Vector{Time},
    tEs::Vector{Time},
)
    ss, hs = residuals_and_waveform_1psr(
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
        tspl,
    )
    spl = CubicHermiteSplineInterpolation([t.t for t in tspl], ss, hs)
    return interp(spl, [t.t for t in tEs])
end

function residuals_1psr_spline(
    mass::Mass,
    n_init::MeanMotion,
    e_init::Eccentricity,
    l0p::InitPhaseParams,
    proj::ProjectionParams,
    dl::Distance,
    dp::Distance,
    α::AzimuthParam,
    z::Redshift,
    terms::Vector{Term},
    tref::Time,
    tEs::Vector{Time},
)
    tspl = spline_time_samples(tEs)
    return residuals_1psr_spline(
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
        tspl,
        tEs,
    )
end
