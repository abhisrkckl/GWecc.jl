export spline_time_samples, residuals_spline, residuals_1psr_spline

struct SimpleHermiteSpline
    xs::Vector{Float64}
    vs::Vector{Float64}
    ts::Vector{Float64}
    tmin::Float64
    tmax::Float64

    function SimpleHermiteSpline(ts, func_and_deriv)
        @assert issorted(ts)
        @assert allunique(ts)
        xs, vs = func_and_deriv(ts)
        new(xs, vs, ts, minimum(ts), maximum(ts))
    end
end

function find_interval(spl::SimpleHermiteSpline, ti::Float64)::Tuple{Int64,Float64,Float64}
    @assert ti >= spl.tmin && ti <= spl.tmax
    for (j, tsplj) in enumerate(spl.ts)
        if ti == tsplj
            return j, tsplj, tsplj
        elseif ti < tsplj
            return j, spl.ts[j-1], tsplj
        end
    end
end

function hermite_spline_basis(u)
    h00 = (1 + 2 * u) * (1 - u)^2
    h10 = u * (1 - u)^2
    h01 = u^2 * (3 - 2 * u)
    h11 = u^2 * (u - 1)
    return h00, h10, h01, h11
end

function interp(spl::SimpleHermiteSpline, ti::Float64)::Float64
    j, t1, t2 = find_interval(spl, ti)

    if t1 == t2
        return spl.xs[j]
    end

    x1 = spl.xs[j-1]
    x2 = spl.xs[j]
    v1 = spl.vs[j-1]
    v2 = spl.vs[j]

    u = (ti - t1) / (t2 - t1)
    h00, h10, h01, h11 = hermite_spline_basis(u)

    dt = t2 - t1

    return h00 * x1 + h10 * dt * v1 + h01 * x2 + h11 * dt * v2
end

function interp(spl::SimpleHermiteSpline, tis::Vector{Float64})::Vector{Float64}
    _interp = ti -> interp(spl, ti)
    return _interp.(tis)
end

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
    terms::Vector{Term},
    tref::Time,
    tspl::Vector{Time},
    tEs::Vector{Time},
)
    res_wav =
        _ts -> residuals_and_waveform(
            mass,
            n_init,
            e_init,
            l0p,
            proj,
            dl,
            dp,
            psrpos,
            gwpos,
            terms,
            tref,
            Time.(_ts),
        )
    spl = SimpleHermiteSpline(extract(tspl), res_wav)
    return interp(spl, extract(tEs))
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
        terms,
        tref,
        tspl,
        tEs,
    )
end

# function residuals_1psr_spline(
#     mass::Mass,
#     n_init::MeanMotion,
#     e_init::Eccentricity,
#     l0p::InitPhaseParams,
#     proj::ProjectionParams,
#     dl::Distance,
#     dp::Distance,
#     α::AzimuthParam,
#     terms::Vector{Term},
#     tref::Time,
#     tspl::Vector{Time},
#     tEs::Vector{Time},
# )
#     res_wav =
#         _ts -> residuals_and_waveform_1psr(
#             mass,
#             n_init,
#             e_init,
#             l0p,
#             proj,
#             dl,
#             dp,
#             α,
#             z,
#             terms,
#             tref,
#             Time.(_ts),
#         )
#     spl = SimpleHermiteSpline(extract(tspl), res_wav)
#     return interp(spl, extract(tEs))
# end

# function residuals_1psr_spline(
#     mass::Mass,
#     n_init::MeanMotion,
#     e_init::Eccentricity,
#     l0p::InitPhaseParams,
#     proj::ProjectionParams,
#     dl::Distance,
#     dp::Distance,
#     α::AzimuthParam,
#     terms::Vector{Term},
#     tref::Time,
#     tEs::Vector{Time},
# )
#     tspl = spline_time_samples(tEs)
#     return residuals_1psr_spline(
#         mass,
#         n_init,
#         e_init,
#         l0p,
#         proj,
#         dl,
#         dp,
#         α,
#         z,
#         terms,
#         tref,
#         tspl,
#         tEs,
#     )
# end
