import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Set plot properties (similar to MATLAB settings)
plt.rcParams['figure.figsize'] = [6, 5]  # Set figure size (width, height)
plt.rcParams['axes.linewidth'] = 0.75    # Axes line width
plt.rcParams['font.size'] = 11           # Font size
plt.rcParams['lines.linewidth'] = 1.5    # Line width
plt.rcParams['lines.markersize'] = 8     # Marker size
plt.rcParams['font.weight'] = 'bold'     # Font weight
plt.rcParams['axes.labelweight'] = 'bold'  # Axis label weight

# Define the physical properties
Tinit = float(input('Enter the Initial Temperature (Kelvin): '))
Tf = float(input('Enter the Final Temperature (Kelvin): '))
Tinf = float(input('Enter the Ambient Temperature (Kelvin): '))
h = float(input('Enter the convection coefficient (h): '))  # convection coefficient (typically 1-20 for air)

# Constants for the iron sphere
ro = 0.1    # 10 cm radius sphere (in meters)
c = 449     # J/kg.K (specific heat capacity of iron)
rho = 7870  # kg/m^3 (density of iron)
k = 2       # W/mK (thermal conductivity of iron)
V = (4/3) * np.pi * ro**3  # Volume of sphere
As = 4 * np.pi * ro**2     # Surface area of sphere
Bi = h * ro / k            # Biot number
Tau = (rho * V * c) / (h * As)  # Time constant

# Time and temperature arrays
ti = 0  # Initial time
Ti = Tinit  # Initial temperature
Tmat = [Ti]  # Array to store temperatures
tmat = [ti]  # Array to store time values
delt = 1     # Time step

if Bi < 0.1:
    # Lumped model case (Bi < 0.1)
    while Tmat[-1] > Tf:
        Tnew = Tinf + (Ti - Tinf) * np.exp(- delt / Tau)
        tnew = tmat[-1] + delt
        Tmat.append(Tnew)
        tmat.append(tnew)
        Ti = Tnew

    # Plot results for lumped model
    plt.plot(tmat, Tmat, 'b-', linewidth=1.5)
    plt.xlabel('Time (s)', fontsize=11, fontweight='bold')
    plt.ylabel('Temperature (K)', fontsize=11, fontweight='bold')
    plt.title('Temperature Variation with Time (Lumped Model)', fontsize=11, fontweight='bold')
    plt.grid(True)
    plt.show()

else:
    # For Bi > 0.1, interpolation is used
    Bimat = np.array([0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0, 40.0, 50.0, 100.0])
    lamat = np.array([0.1730, 0.2445, 0.3450, 0.4217, 0.4860, 0.5423, 0.7593, 0.9208, 1.0528, 1.1656, 1.2644, 1.3525, 1.4320, 1.5044, 1.5708, 2.0288, 2.2889, 2.4556, 2.5704, 2.6537, 2.7165, 2.7654, 2.8044, 2.8363, 2.9857, 3.0372, 3.0632, 3.0788, 3.1102])
    Amat = np.array([1.0030, 1.0060, 1.0120, 1.0179, 1.0239, 1.0298, 1.0592, 1.0880, 1.1164, 1.1441, 1.1713, 1.1978, 1.2236, 1.2488, 1.2732, 1.4793, 1.6227, 1.7202, 1.7870, 1.8383, 1.8673, 1.8920, 1.9106, 1.9249, 1.9781, 1.9898, 1.9942, 1.9962, 1.9990])

    # Interpolate values for Biot number
    A = interp1d(Bimat, Amat, kind='linear', fill_value='extrapolate')(Bi)
    lam = interp1d(Bimat, lamat, kind='linear', fill_value='extrapolate')(Bi)
    alpha = k / (rho * c)  # Thermal diffusivity

    # Minimum time allowed (see Cengel p. 218)
    tmin = int(np.ceil((ro**2) * 0.2 / alpha))
    tmat = np.arange(tmin, 2 * tmin + 1)  # Time values
    Taunew = alpha * tmat / ro**2  # Fourier number

    # Calculate temperature distribution at center of sphere
    Tmatnew = Tinf + (Tinit - Tinf) * A * np.exp(-Taunew * lam**2)

    # Plot results for Bi > 0.1
    plt.plot(tmat, Tmatnew, 'r--', linewidth=1.5)
    plt.xlabel('Time (s)', fontsize=11, fontweight='bold')
    plt.ylabel('Temperature (K)', fontsize=11, fontweight='bold')
    plt.title('Temperature Distribution at Centre of Sphere', fontsize=11, fontweight='bold')
    plt.grid(True)
    
    # Mark region where one-term approximation is valid (Fo > 0.2)
    plt.axvline(tmin, color='k', linestyle=':', label='Fo > 0.2 (One-term Approximation)')
    plt.fill_betweenx([Tf, Tinit], 0, tmin, color='w', alpha=0.5)
    plt.text(tmin / 8, (Tf + Tinit) / 2, 'Fo > 0.2 (One-term Approximation)', fontsize=10)

    # Customize y-axis ticks
    ytickloc = np.round(np.linspace(Tf, Tinit, 10))
    plt.yticks(ytickloc)

    plt.show()
