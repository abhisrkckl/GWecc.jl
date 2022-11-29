export ProjectionParams, SkyLocation
export ScaledTime, Time
export Eccentricity, MeanMotion
export ScaledMeanAnomaly, ScaledPeriastronAngle, Angle, SinCos
export Mass

struct ScaledTime
    τ::Float64
    ScaledTime(τ::Float64) =
        isfinite(τ) && τ >= 0 ? new(τ) : throw(DomainError(τ, "τ<0 encountered."))
end

struct Time
    t::Float64
    Time(t::Float64) = isfinite(t) ? new(t) : throw(DomainError(t, "isnan(t) encountered."))
end

struct Eccentricity
    e::Float64
    Eccentricity(e::Float64) =
        (e > 0 && e < 1) ? new(e) : throw(DomainError(e, "e out of range."))
end

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

struct MeanMotion
    n::Float64
    MeanMotion(n::Float64) =
        (n >= 6.3e-10 && n <= 6.3e-6) ? new(n) : throw(DomainError(n, "n out of range."))
end

struct ScaledMeanAnomaly
    lbar::Float64
    ScaledMeanAnomaly(lbar::Float64) =
        (lbar >= 0 && lbar <= 0.46137) ? new(lbar) :
        throw(DomainError(lbar, "lbar out of range."))
end

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

struct Angle
    θ::Float64
    Angle(θ::Float64) =
        isfinite(θ) ? new(θ) : throw(DomainError(θ, "isnan(θ) encountered."))
end

struct SinCos
    x::Angle
    sinx::Float64
    cosx::Float64

    SinCos(x::Angle) = new(x, sin(x.θ), cos(x.θ))
end

struct ProjectionParams
    ψ::Float64
    cosι::Float64
    γ0::Float64
    γp::Float64

    function ProjectionParams(ψ::Float64, cosι::Float64, γ0::Float64, γp::Float64)
        if !(ψ >= 0 && ψ <= π)
            throw(DomainError(ψ, "ψ out of range."))
        elseif !(cosι >= -1 && cosι < 1)
            throw(DomainError(cosι, "cosι out of range."))
        elseif !(γ0 >= 0 && γ0 <= π)
            throw(DomainError(γ0, "γ0 out of range."))
        elseif !(γp >= 0 && γp <= π)
            throw(DomainError(γp, "γp out of range."))
        else
            return new(ψ, cosι, γ0, γp)
        end
    end
end

struct SkyLocation
    ra::Float64
    dec::Float64

    function SkyLocation(ra::Float64, dec::Float64)
        if !(ra >= 0 && ψ < 2 * π)
            throw(DomainError(ra, "ra out of range."))
        elseif !(dec >= -π / 2 && dec <= π / 2)
            throw(DomainError(dec, "dec out of range."))
        else
            return new(ra, dec)
        end
    end
end
