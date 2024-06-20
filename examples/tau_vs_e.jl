"""Scaled time as a function of eccentricity.
Figure 2 of Susobhanan+ 2020"""

using GWecc
using CairoMakie

println("Running ", PROGRAM_FILE)

es = Eccentricity.(LinRange(0.01, 0.99, 1000))
τs = τ_from_e.(es)

fig = Figure()
ax = Axis(fig[1, 1])
CairoMakie.lines!([τ.τ for τ in τs], [e.e for e in es])
ax.xlabel = "τ"
ax.ylabel = "e"
save("tau_vs_e.pdf", fig)
