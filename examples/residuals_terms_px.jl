"""PTA signals - Earth and Pulsar terms and +/x polarizations"""

using GWecc
using CairoMakie

println("Running ", PROGRAM_FILE)

year = 365.25 * 24 * 3600
parsec = 102927125.0
MSun = 4.92703806e-6

mass = Mass(1e9 * MSun, 0.25)
n_init = MeanMotion(2 * π / (2 * year))
l0p = InitPhaseParams(0.0, 0.0)
tref = Time(5000.0)

proj = ProjectionParams(1e-9, 0.0, 1.0, 0.0, 0.0)

ra_p = 1.5
dec_p = -0.8
ra_gw = 0.5
dec_gw = 0.75
psrpos = SkyLocation(ra_p, dec_p)
gwpos = SkyLocation(ra_gw, dec_gw)

dp = Distance(500 * parsec)

tEs = Time.(LinRange(0, 10 * year, 5000))
tyrs = [t.t for t in tEs] / year

fig = CairoMakie.Figure()

for (idx, e_init) in enumerate(Eccentricity.([0.1, 0.4, 0.7]))

    spEs, sxEs =
        residuals_px(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, EARTH, tref, tEs)
    sEs = residuals(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, [EARTH], tref, tEs)
    spPs, sxPs =
        residuals_px(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, PULSAR, tref, tEs)
    sPs = residuals(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, [PULSAR], tref, tEs)

    ax1 = CairoMakie.Axis(fig[idx, 1])
    CairoMakie.lines!(ax1, tyrs, spEs)
    CairoMakie.lines!(ax1, tyrs, sxEs)
    ax1.ylabel = "s_E(+/x)"
    ax1.xlabel = "t (year)"

    ax2 = CairoMakie.Axis(fig[idx, 2])
    CairoMakie.lines!(ax2, tyrs, sEs)
    ax2.ylabel = "s_E"
    ax2.xlabel = "t (year)"

    ax3 = CairoMakie.Axis(fig[idx, 3])
    CairoMakie.lines!(ax3, tyrs, spPs)
    CairoMakie.lines!(ax3, tyrs, sxPs)
    ax3.ylabel = "s_P(+/x)"
    ax3.xlabel = "t (year)"

    ax4 = CairoMakie.Axis(fig[idx, 4])
    CairoMakie.lines!(ax4, tyrs, sPs)
    ax4.ylabel = "s_P"
    ax4.xlabel = "t (year)"
end

CairoMakie.save("residuals_terms_px.pdf", fig)
