"""
This script contains some external fitting/theoretical functions for
refractive indices
"""
from numpy import sqrt, pi


def AlGaAsIndex(wl, x):
    # from indexAlexey
    e = 1.2395/wl
    e0 = 1.425 + 1.155*x + 0.37 * x**2
    edel = 1.765 + 1.115*x + 0.37 * x**2
    chi = e / e0
    chidel = e / edel
    f1 = (2 - sqrt(1 + chi) - sqrt(1 - chi)) / chi**2
    f2 = (2 - sqrt(1 + chidel) - sqrt(1 - chidel)) / chidel**2
    a1 = 6.3 + 19*x
    b = 9.4 - 10.2*x
    return sqrt(a1*(f1 + 0.5*f2*(e0/edel)**(1.5)) + b)


def SiNxIndex(wl):
    # from Jean Nguyen's Thesis
    C1 = 2.0019336; C2 = 0.15265213; C3 = 4.0495557
    D0=-0.00282; D1=0.003029; D2=-0.0006982; D3=-0.0002839
    D4=0.0001816; D5=-3.948e-005; D6=4.276e-006; D7=-2.314e-007; D8=4.982e-009
    n_SiNx = C1 + C2/wl**2 + C3/wl**4
    k_SiNx = D0+D1*wl+D2*wl**2+D3*wl**3+D4*wl**4+D5*wl**5+D6*wl**6+D7*wl**7+D8*wl**8
    k_SiNx *= 10
    return n_SiNx + 1j*k_SiNx


def SiO2Index(wl):
    # from Jean Nguyen's Thesis
    C1 = 1.41870; C2 = 0.12886725; C3 = 2.7573641e-5
    n_SiO2 = C1 + C2/wl**2 + C3/wl**4
    # this is a 4 peak Lorentzian fit to her data
    y0=-797.4627
    xc1=2.83043; w1=6.083822; A1=10881.9438
    xc2=8.95338; w2=1.38389113; A2=9167.662815
    xc3=12.3845492; w3=3.9792077; A3=12642.72911
    xc4=15.6387213; w4=0.6057751177; A4=3292.325272
    alpha = y0 + 2*A1/pi*w1/(4*(wl-xc1)**2+w1**2) + 2*A2/pi*w2/(4*(wl-xc2)**2+w2**2) \
            + 2*A3/pi*w3/(4*(wl-xc3)**2+w3**2) + 2*A4/pi*w4/(4*(wl-xc4)**2+w4**2)
    k_SiO2 = alpha * wl*1e-4 / (4*pi)
    return n_SiO2 + 1j*k_SiO2
