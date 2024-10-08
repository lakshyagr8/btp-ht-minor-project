import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

# Fin Analysis Function
def fin_analysis(L, r, K):
    h = 35
    Tinf = 20
    Tb = 150
    P = np.pi * 2 * r
    A = (np.pi * r**2) / 4
    Lc = L + ((2 * r) / 4)
    m = np.sqrt((h * P) / (K * A))
    
    xmat = np.arange(0, L / 100, 0.001)
    FImat = Tinf + ((1 / np.exp(m * xmat)) * (Tb - Tinf))
    FAmat = Tinf + ((np.cosh(m * ((L / 100) - xmat)) / np.cosh(m * L / 100)) * (Tb - Tinf))
    FTmat = Tinf + (np.cosh(m * ((Lc / 100) - xmat) / 57.3) / np.cosh((m * Lc / 100) / 57.3) * (Tb - Tinf))
    
    fig, ax = plt.subplots(3, 1, figsize=(6, 8))
    
    ax[0].plot(xmat, FImat, 'r')
    ax[0].set_title('Temperature Distribution for Infinite Fin')
    ax[0].grid(True)

    ax[1].plot(xmat, FAmat, 'b')
    ax[1].set_title('Adiabatic Fin Tip')
    ax[1].grid(True)

    ax[2].plot(xmat, FTmat, 'm')
    ax[2].set_title('Convection at Fin Tip')
    ax[2].grid(True)

    st.pyplot(fig)

# Composite Wall Temperature Profile Function
def composite_wall_temperature_profile(Tin, Tout, L1, K1, L2, K2, L3, K3, A):
    Q = (Tin - Tout) / ((L1 / (K1 * A)) + (L2 / (K2 * A)) + (L3 / (K3 * A)))
    T1 = Tin - (Q * L1 / (K1 * A))
    T2 = T1 - (Q * L2 / (K2 * A))
    T3 = T2 - (Q * L3 / (K3 * A))
    
    st.write(f'Total Heat Transfer through Composite Wall: {Q * 1000:.2f} W')
    
    xmat = np.linspace(0, L1 + L2 + L3, 100)
    ymat = np.piecewise(
        xmat,
        [xmat <= L1, (xmat > L1) & (xmat <= L1 + L2), xmat > L1 + L2],
        [lambda x: Tin - ((Tin - T1) / L1) * x,
         lambda x: T1 - ((T1 - T2) / L2) * (x - L1),
         lambda x: T2 - ((T2 - T3) / L3) * (x - L1 - L2)]
    )

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(xmat, ymat, 'k')
    ax.set_title('Composite Wall Temperature Profile')
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('Temperature (deg C)')
    ax.grid(True)

    st.pyplot(fig)

# Fourier's Law of Heat Conduction Function
def fourier_heat_conduction(X1, X2, T1, T2, A, K):
    delTX = (T2 - T1) / (X2 - X1)  # Temperature Gradient (K/m)
    Q = -K * A * delTX  # Heat transfer rate (W)
    
    st.write(f'THE HEAT TRANSFER RATE BY CONDUCTION IN WATTS IS: {Q:.2f} W')
    
    # Plot
    fig, ax = plt.subplots(figsize=(6, 4))
    y = np.arange(0, 10, 0.1)
    
    if T2 > T1:
        ax.plot([X1] * len(y), y, 'bo', label=f'X1 = {X1}')
        ax.plot([X2] * len(y), y, 'ro', label=f'X2 = {X2}')
    else:
        ax.plot([X1] * len(y), y, 'ro', label=f'X1 = {X1}')
        ax.plot([X2] * len(y), y, 'bo', label=f'X2 = {X2}')
    
    ax.set_xlim([min(X1, X2) - 1, max(X1, X2) + 1])
    ax.set_ylim([0, 10])
    ax.set_xlabel('Position')
    ax.set_ylabel('Y Axis')
    ax.set_title('CONDUCTION HEAT TRANSFER (From Red to Blue)')
    ax.grid(True)
    
    st.pyplot(fig)

# Transient Cooling of Iron Sphere Function
def transient_cooling_of_iron_sphere(Tinit, Tf, Tinf, h):
    ro = 0.1  # radius of the sphere in meters (10 cm)
    c = 449  # J/kg.K for Iron
    rho = 7870  # kg/m^3 for Iron
    k = 2  # W/mK (for Iron)
    V = (4/3) * np.pi * (ro**3)  # volume of sphere
    As = 4 * np.pi * (ro**2)  # surface area of sphere
    Bi = h * ro / k  # Biot number
    Tau = (rho * V * c) / (h * As)  # Time constant

    ti = 0
    Ti = Tinit
    delt = 1  # time step
    Tmat = [Ti]
    tmat = [ti]

    if Bi < 0.1:
        # Lumped system model
        while Tmat[-1] > Tf:
            Tnew = Tinf + (Ti - Tinf) * np.exp(-delt / Tau)
            Ti = Tnew
            tmat.append(tmat[-1] + delt)
            Tmat.append(Tnew)

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(tmat, Tmat, 'b-', linewidth=2)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Temperature (K)')
        ax.set_title('Temperature Variation with Time (Lumped Model)')
        ax.grid(True)
        st.pyplot(fig)

    else:
        st.write("Bi > 0.1: Not suitable for lumped capacitance method.")

# Planck's Law and Wien's Law Function

def planck_wien_law(T):
    # Constants
    C1 = 374200000  # Constant used in Planck law
    C2 = 14390  # Constant used in Planck law
    lambda_max = 2898 / T  # Peak wavelength in microns

    st.write(f'THE PEAK WAVELENGTH (microns) IS {lambda_max:.4f}')

    # Matrices for Power and Wavelength
    ppmat = []
    lpmat = []
    max_power = 0

    # Loop to calculate the data
    for i in np.arange(0.1, 1000, 0.01):
        Lp = i
        pp = C1 / (i ** 5 * (np.exp(C2 / (Lp * T)) - 1))
        ppmat.append(pp)
        lpmat.append(Lp)
        if max_power < pp:
            max_power = pp

    st.write(f'THE PEAK POWER IS {max_power:.4f}')

    # Plotting
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.semilogx(lpmat, ppmat, 'b-', linewidth=1.5)
    ax.axvline(lambda_max, color='r', linestyle='--', label=f'λ_max = {lambda_max:.4f} μm')
    ax.set_xlabel('Wavelength (microns)', fontweight='bold', fontsize=11)
    ax.set_ylabel('Peak Power (W)', fontweight='bold', fontsize=11)
    ax.set_title('Peak Power vs Wavelength', fontweight='bold', fontsize=22)
    ax.grid(True)
    
    # Adjust x-limits based on the scientific notation of lambda_max
    exponent = np.floor(np.log10(lambda_max))
    ax.set_xlim([10 ** (exponent - 1), 10 ** (exponent + 2)])

    ax.legend()
    st.pyplot(fig)



# Streamlit App Interface
st.title('Engineering Analysis Tool')

analysis_type = st.sidebar.selectbox("Choose Analysis", ["Fin Analysis", "Composite Wall Temperature Profile", "Fourier's Law of Heat Conduction", "Transient Cooling of Iron Sphere", "Planck's and Wien's Law"])

if analysis_type == "Fin Analysis":
    L = st.sidebar.number_input('Enter Length (L)', min_value=0.0, value=100.0)
    r = st.sidebar.number_input('Enter Radius (r)', min_value=0.0, value=5.0)
    K = st.sidebar.number_input('Enter Thermal Conductivity (K)', min_value=0.0, value=100.0)
    if st.sidebar.button('Run Fin Analysis'):
        fin_analysis(L, r, K)

elif analysis_type == "Composite Wall Temperature Profile":
    Tin = st.sidebar.number_input('Enter Input Temperature (Tin)', value=28.0)
    Tout = st.sidebar.number_input('Enter Output Temperature (Tout)', value=-30.0)
    L1 = st.sidebar.number_input('Enter Length of Material 1 (L1)', value=8.0)
    K1 = st.sidebar.number_input('Enter Conductivity of Material 1 (K1)', value=0.1)
    L2 = st.sidebar.number_input('Enter Length of Material 2 (L2)', value=20.0)
    K2 = st.sidebar.number_input('Enter Conductivity of Material 2 (K2)', value=0.5)
    L3 = st.sidebar.number_input('Enter Length of Material 3 (L3)', value=15.0)
    K3 = st.sidebar.number_input('Enter Conductivity of Material 3 (K3)', value=3.0)
    A = st.sidebar.number_input('Enter Area (A)', value=1.0)
    if st.sidebar.button('Run Composite Wall Analysis'):
        composite_wall_temperature_profile(Tin, Tout, L1, K1, L2, K2, L3, K3, A)

elif analysis_type == "Fourier's Law of Heat Conduction":
    X1 = st.sidebar.number_input('Enter the value of X1', min_value=0.0, value=0.0)
    X2 = st.sidebar.number_input('Enter the value of X2', min_value=0.0, value=10.0)
    T1 = st.sidebar.number_input('Enter the Temperature of body at X1', value=100.0)
    T2 = st.sidebar.number_input('Enter the Temperature of body at X2', value=50.0)
    A = st.sidebar.number_input('Enter Area of Cross-Section (A)', value=3.0)
    K = st.sidebar.number_input('Enter Thermal Conductivity (K)', value=2.0)
    if st.sidebar.button('Run Fourier Law Analysis'):
        fourier_heat_conduction(X1, X2, T1, T2, A, K)
        
elif analysis_type == "Transient Cooling of Iron Sphere":
    Tinit = st.sidebar.number_input('Enter Initial Temperature (Kelvins)', value=500.0)
    Tf = st.sidebar.number_input('Enter Final Temperature (Kelvins)', value=300.0)
    Tinf = st.sidebar.number_input('Enter Ambient Temperature (Kelvins)', value=293.0)
    h = st.sidebar.number_input('Enter Convective Heat Transfer Coefficient (h)', value=10.0)

    if st.sidebar.button('Run Transient Cooling Analysis'):
        transient_cooling_of_iron_sphere(Tinit, Tf, Tinf, h)

elif analysis_type == "Planck's and Wien's Law":
    T = st.sidebar.number_input('Enter the temperature 1 (in Kelvins)', value=5000.0)
    

    if st.sidebar.button('Run Planck\'s and Wien\'s Law Analysis'):
        planck_wien_law(T)
        
