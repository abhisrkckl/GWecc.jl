include("../src/GWecc.jl")
using .GWecc
using PyCall
using JET

np = pyimport("numpy")

toas = np.linspace(0, 1000000, 100)
pdist = 400.0 # kpc
alpha = 0.3
psi = 1.1
cos_inc = 0.5
log10_M = 9.0
eta = 0.25
log10_F = -8.0
e0 = 0.3
gamma0 = 0.0
gammap = 0.0
l0 = 0.0
lp = 0.0
tref = maximum(toas)
log10_zc = -2.0
psrTerm = false

@report_call eccentric_pta_signal_planck18_1psr(
    toas,
    pdist,
    alpha,
    psi,
    cos_inc,
    log10_M,
    eta,
    log10_F,
    e0,
    gamma0,
    gammap,
    l0,
    lp,
    tref,
    log10_zc,
    psrTerm,
)

@report_opt eccentric_pta_signal_planck18_1psr(
    toas,
    pdist,
    alpha,
    psi,
    cos_inc,
    log10_M,
    eta,
    log10_F,
    e0,
    gamma0,
    gammap,
    l0,
    lp,
    tref,
    log10_zc,
    psrTerm,
)