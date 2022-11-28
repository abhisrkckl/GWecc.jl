export mikkola, kepler, SinCos

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

struct SinCos
    x::Angle
    sinx::Float64
    cosx::Float64
    
    SinCos(x::Angle) = new(x, sin(x.θ), cos(x.θ))
end
