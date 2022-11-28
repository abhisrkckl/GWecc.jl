export mikkola, kepler, true_anomaly_diff, OrbitalPhase

function mikkola(ecc::Eccentricity, ll::Angle)::Angle
    l = ll.θ
    e = ecc.e
    
    if e==0 || l==0 
        return ll
    end

    sgn = sign(l);
    l = sgn*l # l>0
    
    ncycles = floor(l/(2*π));
    l = l - 2*π*ncycles # 0<=l<2*pi

    flag = l > π
    if flag
        l = 2*π - l; # 0<=l<=pi
    end

    alpha = (1-e)/(4*e + 0.5)
    beta = (l/2)/(4*e + 0.5)
    
    z = (beta>0) ? cbrt(beta + sqrt(alpha^3 + beta^2)) : cbrt(beta - sqrt(alpha^3 + beta^2))

    s = (z - alpha/z)
    w = (s - (0.078*s^5)/(1 + e))
    
    E0 = (l + e*(3*w - 4*w^3));
    u = E0
    
    esu  = e*sin(u);
    ecu  = e*cos(u);
    
    fu  = (u - esu - l)
    f1u = (1 - ecu)
    f2u = (esu)
    f3u = (ecu)
    f4u =-(esu)

    u1 = -fu/ f1u
    u2 = -fu/(f1u + f2u*u1/2)
    u3 = -fu/(f1u + f2u*u2/2 + f3u*(u2*u2)/6.0)
    u4 = -fu/(f1u + f2u*u3/2 + f3u*(u3*u3)/6.0 + f4u*(u3*u3*u3)/24.0)
    xi = (E0 + u4)
    
    sol = flag ? (2*π-xi) : xi
    
    u = sgn*(sol + 2*π*ncycles)
    
    return Angle(u)
end

function kepler(ecc::Eccentricity, uu::Angle)::Angle
    u = uu.θ
    e = ecc.e
    return Angle(u - e*sin(u))
end

function true_anomaly_diff(ecc::Eccentricity, scu::SinCos)::Angle
    e = ecc.e
    u = scu.x.θ
    su = scu.sinx
    cu = scu.cosx
    βφ = (e>1e-15) ? (1-sqrt(1-e^2))/e : (e/2 + (e^3)/8 + (e^5)/16)
    v_u = u + 2*atan(βφ*su, 1-βφ*cu)
    v_l = v_u + e*su
    return Angle(v_l)
end

struct OrbitalPhase
    scu::SinCos
    sc2φ::SinCos
    sc2ω::SinCos

    function OrbitalPhase(mass::Mass, n::MeanMotion, et::Eccentricity, l::Angle, γ::Angle)
        scu = SinCos(mikkola(et, l))
        
        eφ = angular_eccentricity(mass, n, et)
        v_l = true_anomaly_diff(eφ, scu).θ
        l = l.θ
        v = v_l + l 
        
        k = advance_of_periastron(mass, n, et).k
        W = (1+k)*v_l

        γ = γ.θ
        φ = γ + l + W
        ω = φ - v

        return new(scu, SinCos(Angle(2*φ)), SinCos(Angle(2*ω)))
    end
end