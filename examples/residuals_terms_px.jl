include("../src/GWecc.jl")
using .GWecc
using PyPlot

year = 365.25 * 24 * 3600
parsec = 102927125.0
MSun = 4.92703806e-6

mass = Mass(1e9 * MSun, 0.25)
n_init = MeanMotion(2 * π / (2 * year))
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

for (idx, e_init) in enumerate(Eccentricity.([0.1, 0.4, 0.7]))

    spEs, sxEs = residuals_px(
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
        EARTH,
        tref,
        tEs,
    )
    sEs = residuals(
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
    spPs, sxPs = residuals_px(
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
        PULSAR,
        tref,
        tEs,
    )
    sPs = residuals(
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

    subplot(3, 4, 4 * (idx - 1) + 1)
    plot(tyrs, spEs)
    plot(tyrs, sxEs)
    ylabel("\$s_{E(+/x)}\$")
    xlabel("t (year)")

    subplot(3, 4, 4 * (idx - 1) + 2)
    plot(tyrs, sEs)
    ylabel("\$s_E\$")
    xlabel("t (year)")

    subplot(3, 4, 4 * (idx - 1) + 3)
    plot(tyrs, spPs)
    plot(tyrs, sxPs)
    ylabel("\$s_{P(+/x)}\$")
    xlabel("t (year)")

    subplot(3, 4, 4 * (idx - 1) + 4)
    plot(tyrs, sPs)
    ylabel("\$s_P\$")
    xlabel("t (year)")
end

show()
