include("../src/GWecc.jl")
using .GWecc
using PyPlot

year = 365.25 * 24 * 3600
parsec = 102927125.0
MSun = 4.92703806e-6

mass = Mass(1e9 * MSun, 0.25)
n_init = MeanMotion(2 * π / (5 * year))
e_init = Eccentricity(0.4)
l_init = Angle(0.0)
tref = Time(5000.0)

proj = ProjectionParams(0.0, 1.0, 0.0, 0.0)

ra_p = 1.5
dec_p = -0.8
ra_gw = 0.5
dec_gw = 0.75
psrpos = SkyLocation(ra_p, dec_p)
gwpos = SkyLocation(ra_gw, dec_gw)

dp = Distance(500 * parsec)
dl = Distance(1e9 * parsec)
z = Redshift(0.1)

tEs = Time.(LinRange(0, 10 * year, 5000))
tyrs = [t.t for t in tEs] / year

𝒜Es = residuals_components_𝒜(
    mass,
    n_init,
    e_init,
    l_init,
    dl,
    dp,
    psrpos,
    gwpos,
    z,
    EARTH,
    tref,
    tEs,
)
𝒜Ps = residuals_components_𝒜(
    mass,
    n_init,
    e_init,
    l_init,
    dl,
    dp,
    psrpos,
    gwpos,
    z,
    PULSAR,
    tref,
    tEs,
)

sEs1 = residuals(
    mass,
    n_init,
    e_init,
    l_init,
    proj,
    dl,
    dp,
    psrpos,
    gwpos,
    z,
    [EARTH],
    tref,
    tEs,
)
sPs1 = residuals(
    mass,
    n_init,
    e_init,
    l_init,
    proj,
    dl,
    dp,
    psrpos,
    gwpos,
    z,
    [PULSAR],
    tref,
    tEs,
)

sEs2 = residuals(
    mass,
    n_init,
    e_init,
    l_init,
    proj,
    dl,
    dp,
    psrpos,
    gwpos,
    z,
    [EARTH],
    tref,
    tEs,
)
sPs2 = residuals(
    mass,
    n_init,
    e_init,
    l_init,
    proj,
    dl,
    dp,
    psrpos,
    gwpos,
    z,
    [PULSAR],
    tref,
    tEs,
)

subplot(221)
for (idx, 𝒜E) in enumerate(𝒜Es)
    plot(tyrs, 𝒜E, label = "\$A_{$idx,E}\$")
end
legend()
ylabel("\$A_{i,E}\$")
xlabel("t (year)")

subplot(222)
plot(tyrs, sEs1)
plot(tyrs, sEs2)
ylabel("\$s_E\$")
xlabel("t (year)")

subplot(223)
for (idx, 𝒜P) in enumerate(𝒜Ps)
    plot(tyrs, 𝒜P, label = "\$A_{$idx,P}\$")
end
legend()
ylabel("\$A_{i,P}\$")
xlabel("t (year)")

subplot(224)
plot(tyrs, sPs1)
plot(tyrs, sPs2)
ylabel("\$s_P\$")
xlabel("t (year)")

show()
