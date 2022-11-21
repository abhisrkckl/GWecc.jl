struct ProjectionParams
    ψ::Float64
    cosι::Float64
    γ0::Float64
    γp::Float64

    function PolarizationAngle(ψ::Float64) 
        if !(ψ>=0 && ψ<=π) 
            throw(DomainError(ψ, "ψ out of range."))
        elseif !(cosι>=-1 && cosι<1)
            throw(DomainError(cosι, "cosι out of range."))
        elseif !(γ0>=0 && γ0<=π) 
            throw(DomainError(γ0, "γ0 out of range."))
        elseif !(γp>=0 && γp<=π) 
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
        if !(ra>=0 && ψ<2*π) 
            throw(DomainError(ra, "ra out of range."))
        elseif !(dec>=-π/2 && dec<=π/2)
            throw(DomainError(dec, "dec out of range."))
        else 
            return new(ra, dec)
        end
    end
end

struct ScaledTime
    τ::Float64
    ScaledTime(τ::Float64) = τ>=0 ? new(τ) : throw(DomainError(τ, "τ<0 encountered."))
end

struct Time
    t::Float64
    Time(t::Float64) = isfinite(t) ? new(t) : throw(DomainError(t, "isnan(t) encountered."))
end

struct Eccentricity
    e::Float64
    Eccentricity(e::Float64) = (e>0 && e<1) ? new(e) : throw(DomainError(e, "e out of range."))
end

struct Mass
    m::Float64
    Mass(m::Float64) = (m>=5e-2 && m<=5e4) ? new(m) : throw(DomainError(m, "m out of range."))
end

struct MeanMotion
    n::Float64
    MeanMotion(n::Float64) = (n>=6.3e-10 && n<=6.3e-6) ? new(n) : throw(DomainError(n, "n out of range."))
end

struct ScaledAngle
    θbar::Float64
    ScaledAngle(θbar::Float64) = (θbar>=0 && θbar<=0.5) ? new(θbar) : throw(DomainError(θbar, "θbar out of range."))
end

struct Angle
    θ::Float64
    ScaledTime(θ::Float64) = θ>=0 ? new(θ) : throw(DomainError(θ, "isnan(θ) encountered."))
end

struct SymmetricMassRatio
    η::Float64
    ScaledAngle(η::Float64) = (η>0 && η<=0.25) ? new(η) : throw(DomainError(η, "η out of range."))
end