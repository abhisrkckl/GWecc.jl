export residual_spline, residual_1psr_spline

using CubicHermiteSpline

function spline_time_samples(tEs::Vector{Time})
    day = 24 * 3600
    tEs_day = tEs / day
    return unique(round.(tEs_day)) * day
end

function residual_spline(
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
    hs, ss = residuals_and_waveform(
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

function residual_1psr_spline(
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
    hs, ss = residuals_and_waveform_1psr(
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


