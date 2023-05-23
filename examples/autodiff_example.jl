using GWecc
using Zygote
using PyPlot

year = 365.25 * 24 * 3600
toas = LinRange(0, 10 * year, 5)

theta = π / 3
phi = π / 4
psrdist = 1.0
cos_gwtheta = 0.3
gwphi = π / 5
psi = 0.0
cos_inc = 0.5
log10_M = 8.5
eta = 0.2
log10_F = -8.0
e0 = 0.3
gamma0 = gammap = 0.0
l0 = lp = 0.0
tref = maximum(toas)
log10_A = -8.0

psrTerm = false
spline = false

eccentric_pta_signal_for_ad(toas, tref, psrTerm, spline) =
    (
        theta,
        phi,
        psrdist,
        cos_gwtheta,
        gwphi,
        psi,
        cos_inc,
        log10_M,
        eta,
        log10_F,
        e0,
        gamma0,
        l0,
        log10_A,
    ) -> eccentric_pta_signal(
        toas,
        theta,
        phi,
        psrdist,
        cos_gwtheta,
        gwphi,
        psi,
        cos_inc,
        log10_M,
        eta,
        log10_F,
        e0,
        gamma0,
        gamma0,
        l0,
        l0,
        tref,
        log10_A,
        psrTerm,
        spline,
    )

gwe = eccentric_pta_signal_for_ad(toas, tref, psrTerm, spline)
Rs = gwe(
    theta,
    phi,
    psrdist,
    cos_gwtheta,
    gwphi,
    psi,
    cos_inc,
    log10_M,
    eta,
    log10_F,
    e0,
    gamma0,
    l0,
    log10_A,
)

PyPlot.plot(toas, Rs)
PyPlot.show()

J = jacobian(
    gwe,
    theta,
    phi,
    psrdist,
    cos_gwtheta,
    gwphi,
    psi,
    cos_inc,
    log10_M,
    eta,
    log10_F,
    e0,
    gamma0,
    l0,
    log10_A,
)
