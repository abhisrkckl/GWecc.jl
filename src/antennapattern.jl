export sky_direction_uvec, gw_polarization_tensors, AntennaPattern, pulsar_term_delay

using LinearAlgebra

"Unit vector pointing to a sky location."
function sky_direction_uvec(pos::SkyLocation)
    λ = pos.ra
    β = pos.dec

    n1 = cos(λ) * cos(β)
    n2 = sin(λ) * cos(β)
    n3 = sin(β)

    return [n1, n2, n3]
end

"+/x GW polarization tensors corresponding to a sky location"
function gw_polarization_tensors(pos::SkyLocation)
    λ = pos.ra
    β = pos.dec
    e11p = (sin(λ)^2) - (cos(λ)^2) * (sin(β)^2)
    e21p = -sin(λ) * cos(λ) * ((sin(β)^2) + 1)
    e31p = cos(λ) * sin(β) * cos(β)
    e12p = -sin(λ) * cos(λ) * ((sin(β)^2) + 1)
    e22p = (cos(λ)^2) - (sin(λ)^2) * (sin(β)^2)
    e32p = sin(λ) * sin(β) * cos(β)
    e13p = cos(λ) * sin(β) * cos(β)
    e23p = sin(λ) * sin(β) * cos(β)
    e33p = -(cos(β)^2)

    e11c = sin(2 * λ) * sin(β)
    e21c = -cos(2 * λ) * sin(β)
    e31c = -sin(λ) * cos(β)
    e12c = -cos(2 * λ) * sin(β)
    e22c = -sin(2 * λ) * sin(β)
    e32c = cos(λ) * cos(β)
    e13c = -sin(λ) * cos(β)
    e23c = cos(λ) * cos(β)
    e33c = 0

    ep = [e11p e12p e13p; e21p e22p e23p; e31p e32p e33p]
    ec = [e11c e12c e13c; e21c e22c e23c; e31c e32c e33c]

    return ep, ec
end

"""Antenna pattern functions Fp and Fx and the cos-angle 
between the pulsar and source locations."""
struct AntennaPattern
    cosµ::Float64
    Fp::Float64
    Fx::Float64

    function AntennaPattern(psrpos::SkyLocation, gwpos::SkyLocation)
        # cosµ = cos(β) * cos(βp) * cos(λ - λp) + sin(β) * sin(βp)

        phat = sky_direction_uvec(psrpos)
        nhat = sky_direction_uvec(gwpos)
        ep, ex = gw_polarization_tensors(gwpos)

        cosµ = dot(phat, nhat)
        if isapprox(cosµ, 1.0, atol = 1e-6)
            throw(DomainError(cosµ, "cosµ ≈ 1 encountered."))
        end

        Fp = dot(phat, ep, phat) * 0.5 / (1 - cosµ)
        Fx = dot(phat, ex, phat) * 0.5 / (1 - cosµ)

        return new(cosµ, Fp, Fx)
    end
end

"Time delay between the earth and the pulsar terms."
function pulsar_term_delay(ap::AntennaPattern, psrdist::Distance)::Time
    dp = psrdist.D
    cosµ = ap.cosµ
    return Time(-dp * (1 - cosµ))
end
