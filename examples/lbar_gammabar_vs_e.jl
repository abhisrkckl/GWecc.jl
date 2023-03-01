"""Scaled mean anomaly and scaled periastron angle as functions of eccentricity.
Figure 3 of Susobhanan+ 2020"""

using GWecc
using PyPlot

println("Running ", PROGRAM_FILE)

es = Eccentricity.(LinRange(0.01, 0.99, 1000))
lbars = lbar_from_e.(es)
γbars = γbar_from_e.(es)

plot([e.e for e in es], [lbar.lbar for lbar in lbars])
plot([e.e for e in es], γbars)
ylabel("\$\\bar{l}\$, \$\\bar{\\gamma}\$")
xlabel("\$e\$")
show()
