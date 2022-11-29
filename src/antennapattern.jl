export AntennaPattern

struct AntennaPattern
    cosµ::Float64
    Fp::Float64
    Fx::Float64

    function AntennaPattern(psrpos::SkyLocation, gwpos::SkyLocation)
        λp = psrpos.ra
        βp = psrpos.dec
        λ = gwpos.ra
        β = gwpos.dec

        cosµ = cos(β) * cos(βp) * cos(λ - λp) + sin(β) * sin(βp)

        if isapprox(cosµ, 1.0, atol = 1e-6)
            throw(DomainError(cosµ, "cosµ ≈ 1 encountered."))
        end

        n1 = cos(λp) * cos(βp)
        n2 = sin(λp) * cos(βp)
        n3 = sin(βp)

        e11p = (sin(λ)^2) - (cos(λ)^2) * (sin(β)^2)
        e21p = -sin(λ) * cos(λ) * ((sin(β)^2) + 1)
        e31p = cos(λ) * sin(β) * cos(β)
        e12p = -sin(λ) * cos(λ) * ((sin(β)^2) + 1)
        e22p = (cos(λ)^2) - (sin(λ)^2) * (sin(β)^2)
        e32p = sin(λ) * sin(β) * cos(β)
        e13p = cos(λ) * sin(β) * cos(β)
        e23p = sin(λ) * sin(β) * cos(β)
        e33p = -(cos(β)^2)

        Fp =
            (
                n1 * (n1 * e11p + n2 * e12p + n3 * e13p) +
                n2 * (n1 * e21p + n2 * e22p + n3 * e23p) +
                n3 * (n1 * e31p + n2 * e32p + n3 * e33p)
            ) * 0.5 / (1 - cosµ)

        e11c = sin(2 * λ) * sin(β)
        e21c = -cos(2 * λ) * sin(β)
        e31c = -sin(λ) * cos(β)
        e12c = -cos(2 * λ) * sin(β)
        e22c = -sin(2 * λ) * sin(β)
        e32c = cos(λ) * cos(β)
        e13c = -sin(λ) * cos(β)
        e23c = cos(λ) * cos(β)
        e33c = 0

        Fx =
            (
                n1 * (n1 * e11c + n2 * e12c + n3 * e13c) +
                n2 * (n1 * e21c + n2 * e22c + n3 * e23c) +
                n3 * (n1 * e31c + n2 * e32c + n3 * e33c)
            ) * 0.5 / (1 - cosµ)

        return new(cosµ, Fp, Fx)
    end
end
