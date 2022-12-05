"""Comparison of eccentric waveforms and residuals (+/x polarizations) for different eccentricities.
"""

include("../src/GWecc.jl")
using .GWecc
using PyPlot

year = 365.25 * 24 * 3600
parsec = 102927125.0
MSun = 4.92703806e-6

mass = Mass(1e9 * MSun, 0.25)
n_init = MeanMotion(2 * π / (2 * year))
l_init = Angle(0.0)

dl = Distance(1e9 * parsec)
proj = ProjectionParams(0.0, 1.0, 0.0, 0.0)

ts = Time.(LinRange(0, 10 * year, 5000))

for (idx, e_init) in enumerate(Eccentricity.([0.1, 0.4, 0.8]))
    coeffs = EvolvCoeffs(mass, n_init, e_init)
    
    hpxs = [waveform_px(mass, coeffs, l_init, proj, dl, false, dt) for dt in ts]
    hps = [hpx[1] for hpx in hpxs]
    hxs = [hpx[2] for hpx in hpxs]
    
    subplot(320 + 2*idx - 1)
    plot([t.t for t in ts]/year, hps)
    plot([t.t for t in ts]/year, hxs)
    ylabel("\$h_{+,\\times}\$")
    xlabel("t (year)")

    spxs = [residual_px(mass, coeffs, l_init, proj, dl, false, dt) for dt in ts]
    sps = [spx[1] for spx in spxs]
    sxs = [spx[2] for spx in spxs]

    subplot(320 + 2*idx)
    plot([t.t for t in ts]/year, sps)
    plot([t.t for t in ts]/year, sxs)
    ylabel("\$s_{+,\\times}\$")
    xlabel("t (year)")
end

show()