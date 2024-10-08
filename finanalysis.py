import numpy as np
import matplotlib.pyplot as plt

def fin_analysis(L, r, K):
    h = 35
    Tinf = 20
    Tb = 150
    P = np.pi * 2 * r
    A = (np.pi * r**2) / 4
    Lc = L + ((2 * r) / 4)
    m = np.sqrt((h * P) / (K * A))
    M = np.sqrt(h * P * K * A)
    
    xmat = np.arange(0, L / 100, 0.001)
    FImat = Tinf + ((1 / np.exp(m * xmat)) * (Tb - Tinf))
    FAmat = Tinf + ((np.cosh(m * ((L / 100) - xmat)) / np.cosh(m * L / 100)) * (Tb - Tinf))
    FTmat = Tinf + (np.cosh(m * ((Lc / 100) - xmat) / 57.3) / np.cosh((m * Lc / 100) / 57.3) * (Tb - Tinf))
    
    plt.figure(figsize=(6, 5))

    plt.subplot(3, 1, 1)
    plt.plot(xmat, FImat, 'r', linewidth=1.5)
    plt.title('TEMPERATURE DISTRIBUTION FOR INFINITE FIN')
    plt.xlabel('x (mm)')
    plt.ylabel('Temperature T (deg C)')
    plt.grid(True)

    plt.subplot(3, 1, 2)
    plt.plot(xmat, FAmat, 'b', linewidth=1.5)
    plt.title('TEMPERATURE DISTRIBUTION WITH ADIABATIC FIN TIP')
    plt.xlabel('x (mm)')
    plt.ylabel('Temperature T (deg C)')
    plt.grid(True)

    plt.subplot(3, 1, 3)
    plt.plot(xmat, FTmat, 'm', linewidth=1.5)
    plt.title('TEMPERATURE DISTRIBUTION WITH CONVECTION AT FIN TIP')
    plt.xlabel('x (mm)')
    plt.ylabel('Temperature T (deg C)')
    plt.grid(True)

    plt.tight_layout()
    plt.show()
