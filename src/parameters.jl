export ProjectionParams, SkyLocation, InitPhaseParams
export ScaledTime, Time, Distance, Redshift, unredshifted_time_difference
export Eccentricity, MeanMotion
export ScaledMeanAnomaly, ScaledPeriastronAngle, Angle, SinCos
export Mass
export Term, EARTH, PULSAR

import Base.+, Base.-, Base.*

"Dimensionless scaled time (Defined in Susobhanan+ 2020, Section IIC)"
struct ScaledTime
    τ::Float64
    ScaledTime(τ::Float64) =
        isfinite(τ) && τ >= 0 ? new(τ) : throw(DomainError(τ, "τ<0 encountered."))
end

"Time in seconds"
struct Time
    t::Float64
    Time(t::Float64) = isfinite(t) ? new(t) : throw(DomainError(t, "isnan(t) encountered."))
end
t1::Time + t2::Time = Time(t1.t + t2.t)
t1::Time - t2::Time = Time(t1.t - t2.t)
-t1::Time = Time(-t1.t)
a::Number * t1::Time = Time(a * t1.t)
t1::Time * a::Number = a * t1

"Eccentricity. Must lie in [0,1)."
struct Eccentricity
    e::Float64
    Eccentricity(e::Float64) =
        (e > 0 && e < 1) ? new(e) : throw(DomainError(e, "e out of range."))
end

"""Total mass, symmetric mass ratio and chirp mass. 
Masses are represented in geometric units (s).
Symmetric mass ratio must lie in (0, 0.25]."""
struct Mass
    m::Float64
    η::Float64
    Mch::Float64
    function Mass(m::Float64, η::Float64)
        if !(m >= 5e-2 && m <= 5e4)
            throw(DomainError(m, "m out of range."))
        elseif !(η > 0 && η <= 0.25)
            throw(DomainError(η, "η out of range."))
        else
            Mc = m * η^(3 / 5)
            new(m, η, Mc)
        end
    end
end

"Mean motion in rad/s"
struct MeanMotion
    n::Float64
    MeanMotion(n::Float64) =
        (n >= 1e-11 && n <= 5e-6) ? new(n) : throw(DomainError(n, "n out of range."))
end

"Dimensionless scaled mean anomaly."
struct ScaledMeanAnomaly
    lbar::Float64
    ScaledMeanAnomaly(lbar::Float64) =
        (lbar >= 0 && lbar <= 0.46137) ? new(lbar) :
        throw(DomainError(lbar, "lbar out of range."))
end

"Dimensionless scaled periastron angles at the 1PN, 2PN and 3PN orders."
struct ScaledPeriastronAngle
    γbar1::Float64
    γbar2::Float64
    γbar3::Float64
    function ScaledPeriastronAngle(γbar1::Float64, γbar2::Float64, γbar3::Float64)
        if !(γbar1 > 0 && γbar1 < 0.084876)
            throw(DomainError(γbar1, "γbar1 out of range."))
        elseif !(γbar2 > 0 && γbar2 < 2.4914)
            throw(DomainError(γbar2, "γbar2 out of range."))
        elseif !(isfinite(γbar3) && γbar3 < -43.2093)
            throw(DomainError(γbar3, "γbar3 out of range."))
        else
            new(γbar1, γbar2, γbar3)
        end
    end
end

"Angle in rad"
struct Angle
    θ::Float64
    Angle(θ::Float64) =
        isfinite(θ) ? new(θ) : throw(DomainError(θ, "isnan(θ) encountered."))
end

"Initial mean anomalies for the earth and the pulsar terms."
struct InitPhaseParams
    l0::Angle
    lp::Angle
    function InitPhaseParams(l0::Float64, lp::Float64)
        if l0 < 0 || l0 >= 2 * π
            throw(DomainError(l0, "l0 out of range."))
        elseif lp < 0 || lp >= 2 * π
            throw(DomainError(lp, "lp out of range."))
        else
            return new(Angle(l0), Angle(lp))
        end
    end
end

"sin and cos of an angle."
struct SinCos
    x::Angle
    sinx::Float64
    cosx::Float64

    SinCos(x::Angle) = new(x, sin(x.θ), cos(x.θ))
end

"""Projection parameters including polarization angle, 
inclination, and the initial periastron angles of the 
earth and the pulsar terms."""
struct ProjectionParams
    sc2ψ::SinCos
    cosι::Float64
    γ0::Float64
    γp::Float64

    function ProjectionParams(ψ::Float64, cosι::Float64, γ0::Float64, γp::Float64)
        if !(ψ >= 0 && ψ < π)
            throw(DomainError(ψ, "ψ out of range."))
        elseif !(cosι >= -1 && cosι <= 1)
            throw(DomainError(cosι, "cosι out of range."))
        elseif !(γ0 >= 0 && γ0 <= π)
            throw(DomainError(γ0, "γ0 out of range."))
        elseif !(γp >= 0 && γp <= π)
            throw(DomainError(γp, "γp out of range."))
        else
            return new(SinCos(Angle(2 * ψ)), cosι, γ0, γp)
        end
    end
end

"Sky location represented by RA and DEC in rad"
struct SkyLocation
    ra::Float64
    dec::Float64

    function SkyLocation(ra::Float64, dec::Float64)
        if !(ra >= 0 && ra < 2 * π)
            throw(DomainError(ra, "ra out of range."))
        elseif !(dec >= -π / 2 && dec <= π / 2)
            throw(DomainError(dec, "dec out of range."))
        else
            return new(ra, dec)
        end
    end
end

"Distance in geometric units (s)"
struct Distance
    D::Float64
    Distance(D::Float64) =
        D > 0 ? new(D) : throw(DomainError(D, "distance should be positive."))
end

"Cosmological redshift. Must be positive."
struct Redshift
    z::Float64
    Redshift(z::Float64) =
        z >= 0 ? new(z) : throw(DomainError(z, "redshift should be positive."))
end

"Apply redshift to a time difference"
unredshifted_time_difference(t::Time, tref::Time, z::Redshift)::Time = (t - tref) / (1 + z.z)

@enum Term EARTH PULSAR
