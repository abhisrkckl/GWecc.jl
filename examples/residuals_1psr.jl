import PyPlot
using GWecc

println("Running ", PROGRAM_FILE)

mass = Mass(5000.0, 0.1)
n_init = MeanMotion(1e-8)
e_init = Eccentricity(0.5)
l_init = Angle(0.1)
γ_init = Angle(1.25)
dl = Distance(1e16)

ra_p = 0.0
dec_p = pi / 6
ra_gw = pi / 2
dec_gw = 0.0
psrpos = SkyLocation(ra_p, dec_p)
gwpos = SkyLocation(ra_gw, dec_gw)
dp = Distance(1e13)

ψ = 0.0
cosι = 0.0
γ0 = 0.5
γp = γ0

S0 = 1e-9

proj = ProjectionParams(S0, ψ, cosι, γ0)

l0p = InitPhaseParams(l_init.θ)

year = 365.25 * 24 * 3600
tref = Time(60 * year)
tEs = time_range(Time(0.0), tref, 1000)

ap = AntennaPattern(psrpos, gwpos)
Δp = pulsar_term_delay(ap, dp)

proj1 = ProjectionParams1psr(proj, ap)

Rs = residuals(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, [EARTH], tref, tEs)
Rs1 = residuals_1psr(mass, n_init, e_init, l_init, proj1, Δp, [EARTH], tref, tEs)

PyPlot.plot(extract(tEs)/year, Rs*1e9, label="Full waveform")
PyPlot.plot(extract(tEs)/year, Rs1*1e9, label="Single pulsar")
PyPlot.legend()
PyPlot.xlabel("t (year)")
PyPlot.ylabel("R(t) (ns)")
PyPlot.show()

