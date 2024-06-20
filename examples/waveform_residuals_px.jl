"""Comparison of eccentric waveforms and residuals (+/x polarizations) for different eccentricities."""

using GWecc
using CairoMakie

println("Running ", PROGRAM_FILE)

year = 365.25 * 24 * 3600
parsec = 102927125.0
MSun = 4.92703806e-6

mass = Mass(1e9 * MSun, 0.25)
n_init = MeanMotion(2 * π / (2 * year))
l0p = InitPhaseParams(0.0, 0.0)

proj = ProjectionParams(1e-9, 0.0, 1.0, 0.0, 0.0)

ts = Time.(LinRange(0, 10 * year, 5000))

fig = CairoMakie.Figure()

for (idx, e_init) in enumerate(Eccentricity.([0.1, 0.4, 0.8]))
    coeffs = EvolvCoeffs(mass, n_init, e_init)

    hpxs = [waveform_px(mass, coeffs, l0p, proj, false, dt) for dt in ts]
    hps = [hpx[1] for hpx in hpxs]
    hxs = [hpx[2] for hpx in hpxs]

    ax1 = CairoMakie.Axis(fig[idx, 1])
    CairoMakie.lines!(ax1, [t.t for t in ts] / year, hps)
    CairoMakie.lines!(ax1, [t.t for t in ts] / year, hxs)
    ax1.ylabel = "h+, hx"
    ax1.xlabel = "t (year)"

    spxs = [residual_px(mass, coeffs, l0p, proj, false, dt) for dt in ts]
    sps = [spx[1] for spx in spxs]
    sxs = [spx[2] for spx in spxs]

    ax2 = CairoMakie.Axis(fig[idx, 2])
    CairoMakie.lines!(ax2, [t.t for t in ts] / year, sps)
    CairoMakie.lines!(ax2, [t.t for t in ts] / year, sxs)
    ax2.ylabel = "s+, sx"
    ax2.xlabel = "t (year)"
end

save("waveform_residuals_px.pdf", fig)
