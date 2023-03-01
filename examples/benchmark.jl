using GWecc
using BenchmarkTools

mass = Mass(5000.0, 0.1)
n_init = MeanMotion(1e-8)
e_init = Eccentricity(0.1)
l_init = Angle(0.1)
γ_init = Angle(1.25)

ra_p = 1.5
dec_p = -0.8
ra_gw = 0.5
dec_gw = 0.75
psrpos = SkyLocation(ra_p, dec_p)
gwpos = SkyLocation(ra_gw, dec_gw)
dp = Distance(1e13)
ap = AntennaPattern(psrpos, gwpos)

ψ = 1.1
cosι = 0.52
γ0 = γ_init.θ
γp = γ0 + 0.2
S0 = 1e-9
proj = ProjectionParams(S0, ψ, cosι, γ0, γp)

l0p = InitPhaseParams(l_init.θ, l_init.θ)

year = 365.25 * 24 * 3600
_tref = 10 * year
tref = Time(_tref)
_tEs = LinRange(0, _tref, 100)
_tEs = reduce(vcat, [t .+ LinRange(0, 2, 16) for t in _tEs])
tEs = Time.(_tEs)

term = EARTH

hs = waveform(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, [term], tref, tEs)
rs = residuals(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, [term], tref, tEs)
rs, hs = residuals_and_waveform(
    mass,
    n_init,
    e_init,
    l0p,
    proj,
    dp,
    psrpos,
    gwpos,
    [term],
    tref,
    tEs,
)
rs_spl =
    residuals_spline(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, [term], tref, tEs)

print("waveform :: ")
@btime waveform(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, [term], tref, tEs)

print("residuals :: ")
@btime residuals(mass, n_init, e_init, l0p, proj, dp, psrpos, gwpos, [term], tref, tEs)

print("residuals_and_waveform :: ")
@btime residuals_and_waveform(
    mass,
    n_init,
    e_init,
    l0p,
    proj,
    dp,
    psrpos,
    gwpos,
    [term],
    tref,
    tEs,
)

print("residuals_spline :: ")
@btime residuals_spline(
    mass,
    n_init,
    e_init,
    l0p,
    proj,
    dp,
    psrpos,
    gwpos,
    [term],
    tref,
    tEs,
)
