"""Scaled mean anomaly and scaled periastron angle as functions of eccentricity.
Figure 3 of Susobhanan+ 2020"""

using GWecc
import GLMakie

println("Running ", PROGRAM_FILE)

es = Eccentricity.(LinRange(0.01, 0.99, 1000))
lbars = lbar_from_e.(es)
γbars = γbar_from_e.(es)

fig, ax, plt = GLMakie.lines([e.e for e in es], [lbar.lbar for lbar in lbars])
GLMakie.lines!([e.e for e in es], γbars)
ax.ylabel = "lbar, gammabar"
ax.xlabel = "e"
GLMakie.current_figure()
GLMakie.save("lbar_gammabar_vs_e.png", fig)