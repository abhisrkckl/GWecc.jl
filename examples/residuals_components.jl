"""PTA signal components for fast likelihood"""

using GWecc
using CairoMakie

println("Running ", PROGRAM_FILE)

year = 365.25 * 24 * 3600
parsec = 102927125.0
MSun = 4.92703806e-6

mass = Mass(1e9 * MSun, 0.25)
n_init = MeanMotion(2 * π / (5 * year))
e_init = Eccentricity(0.4)
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

𝒜Es = residuals_components_𝒜(mass, n_init, e_init, l0p, dp, psrpos, gwpos, EARTH, tref, tEs)
𝒜Ps =
    residuals_components_𝒜(mass, n_init, e_init, l0p, dp, psrpos, gwpos, PULSAR, tref, tEs)

sEs1 = residuals(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, [EARTH], tref, tEs)
sPs1 = residuals(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, [PULSAR], tref, tEs)

sEs2 = residuals(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, [EARTH], tref, tEs)
sPs2 = residuals(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, [PULSAR], tref, tEs)

fig = CairoMakie.Figure()
ax1 = CairoMakie.Axis(fig[1,1])
for (idx, 𝒜E) in enumerate(𝒜Es)
    CairoMakie.lines!(ax1, tyrs, 𝒜E, label = "A_$idx[E]")
end
CairoMakie.axislegend(ax1, framevisible=false)
ax1.ylabel = "A_i[E]"
ax1.xlabel = "t (year)"

ax2 = CairoMakie.Axis(fig[1,2])
CairoMakie.lines!(ax2, tyrs, sEs1)
CairoMakie.lines!(ax2, tyrs, sEs2)
ax2.ylabel = "s_E"
ax2.xlabel = "t (year)"

ax3 = CairoMakie.Axis(fig[2,1])
for (idx, 𝒜P) in enumerate(𝒜Ps)
    CairoMakie.lines!(ax3, tyrs, 𝒜P, label = "A_$idx[P]")
end
CairoMakie.axislegend(ax3, framevisible=false)
ax3.ylabel = "A_i[P]"
ax3.xlabel = "t (year)"

ax4 = CairoMakie.Axis(fig[2,2])
CairoMakie.lines!(ax4, tyrs, sPs1)
CairoMakie.lines!(tyrs, sPs2)
ax4.ylabel = "s_P"
ax4.xlabel = "t (year)"

CairoMakie.save("residuals_components.pdf", fig)
