import numpy as np

def f(t, τ, γ, vy, r0, wpar, wper, α):
    return ((4 - γ**2) / γ**2 * t/τ - 4*np.sin(γ*vy*t/r0)/γ**3 * r0/(vy*τ)
            - 2*(1-np.cos(γ*vy*t/r0)) / γ**2 * wpar/wper * r0/(vy*τ) * np.sin(α))

def g(t, τ, γ, vy, r0, wpar, wper, α, b):
    w = np.sqrt(wpar**2 + wper**2)
    return 2*(1-np.cos(γ*vy*t/r0)) / γ**2 * b*w**2*np.cos(α) / (r0*wper**2) * r0/(vy*τ)

def B2(b, rs, r0, wpar, wper):
    w = np.sqrt(wpar**2 + wper**2)
    return (b**2 + rs**2) / r0**2 * w**2 / wper**2

def rho(t, f, B2, psi0, g):
    return (1 + (f*B2 - f*psi0**2 + 2*g*psi0) / (psi0**2 + B2)**2) ** (-1.)
