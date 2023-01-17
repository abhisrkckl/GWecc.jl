"""Scaled time as a function of eccentricity.
Figure 2 of Susobhanan+ 2020"""

using GWecc
using PyPlot

es = Eccentricity.(LinRange(0.01, 0.99, 1000))
τs = τ_from_e.(es)

plot([τ.τ for τ in τs], [e.e for e in es])
xlabel("τ")
ylabel("e")
show()
