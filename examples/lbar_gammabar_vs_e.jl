# Figure 3 of Susobhanan+ 2020

include("../src/GWecc.jl")
using .GWecc
using PyPlot

es = Eccentricity.(LinRange(0.01, 0.99, 1000))
lbars = lbar_from_e.(es)
γbars = γbar_from_e.(es)

plot([e.e for e in es], [lbar.lbar for lbar in lbars])
plot([e.e for e in es], γbars)
ylabel("\$\\bar{l}\$, \$\\bar{\\gamma}\$")
xlabel("\$e\$")
show()