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

def transient_cooling_with_material_selection(Tinit, Tf, Tinf, h, material):
    materials = {
    'Iron': {'c': 449, 'rho': 7870, 'k': 80},  # J/kg.K, kg/m^3, W/mK
    'Copper': {'c': 385, 'rho': 8960, 'k': 401},
    'Aluminum': {'c': 900, 'rho': 2700, 'k': 237},
    'Silver': {'c': 235, 'rho': 10490, 'k': 429},
    'Gold': {'c': 129, 'rho': 19320, 'k': 318},
    'Brass': {'c': 380, 'rho': 8530, 'k': 109},
    'Nickel': {'c': 444, 'rho': 8908, 'k': 91},
    'Lead': {'c': 128, 'rho': 11340, 'k': 35},
    'Zinc': {'c': 388, 'rho': 7130, 'k': 116},
    'Titanium': {'c': 523, 'rho': 4506, 'k': 22},
    'Tungsten': {'c': 134, 'rho': 19300, 'k': 173},
    'Platinum': {'c': 133, 'rho': 21450, 'k': 72},
}
    c = materials[material]['c']
    rho = materials[material]['rho']
    k = materials[material]['k']
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
        while Tmat[-1] > Tf:
            Tnew = Tinf + (Ti - Tinf) * np.exp(-delt / Tau)
            Ti = Tnew
            tmat.append(tmat[-1] + delt)
            Tmat.append(Tnew)

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(tmat, Tmat, 'b-', linewidth=2)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Temperature (K)')
        ax.set_title(f'Temperature Variation with Time (Lumped Model, {material})')
        ax.grid(True)
        st.pyplot(fig)

    else:
        st.write("Bi > 0.1: Not suitable for lumped capacitance method.")

# Planck's Law and Wien's Law with Multiple Graphs
def planck_wien_law_multi(T_values):
    C1 = 374200000  # Constant used in Planck law
    C2 = 14390  # Constant used in Planck law
    colors = ['b', 'g', 'r']  # Colors for the graphs

    fig, ax = plt.subplots(figsize=(6, 4))

    for i, T in enumerate(T_values):
        lambda_max = 2898 / T  # Peak wavelength in microns
        st.write(f'THE PEAK WAVELENGTH FOR T={T} K IS {lambda_max:.4f} MICRONS')

        ppmat = []
        lpmat = []

        for Lp in np.arange(0.1, 1000, 0.01):
            pp = C1 / (Lp ** 5 * (np.exp(C2 / (Lp * T)) - 1))
            ppmat.append(pp)
            lpmat.append(Lp)

        ax.semilogx(lpmat, ppmat, colors[i % len(colors)], linewidth=1.5, label=f'T = {T} K')
        ax.axvline(lambda_max, color=colors[i % len(colors)], linestyle='--', label=f'λ_max for {T} K')

    ax.set_xlabel('Wavelength (microns)', fontweight='bold', fontsize=11)
    ax.set_ylabel('Power (W)', fontweight='bold', fontsize=11)
    ax.set_title('Power vs Wavelength', fontweight='bold', fontsize=16)
    ax.grid(True)
    ax.legend()
    st.pyplot(fig)

# 3-Body Enclosure Radiation Problem
def three_body_radiation(T1, T2, T3, e1, e2, e3, A1, A2, A3):
    # Stefan-Boltzmann constant
    sigma = 5.67e-8  # W/m².K⁴

    # Net radiation exchange (assuming view factors Fij = 1 between all surfaces)
    Q12 = sigma * A1 * (e1 / (1 - e1)) * (T1**4 - T2**4)
    Q13 = sigma * A1 * (e1 / (1 - e1)) * (T1**4 - T3**4)
    Q23 = sigma * A2 * (e2 / (1 - e2)) * (T2**4 - T3**4)

    st.write(f'Net Heat Exchange Between Body 1 and 2: {Q12:.2f} W')
    st.write(f'Net Heat Exchange Between Body 1 and 3: {Q13:.2f} W')
    st.write(f'Net Heat Exchange Between Body 2 and 3: {Q23:.2f} W')

# Temperature and Velocity Profile for Flow Over a Flat Plate
def flat_plate_profiles(Uinf, Tw, Tinf, x, Pr):
    nu = 1.5e-5  # Kinematic viscosity (m^2/s) for air at room temp
    alpha = nu / Pr  # Thermal diffusivity

    # Boundary layer thickness and thermal boundary layer thickness
    delta = 5 * np.sqrt(nu * x / Uinf)
    delta_t = 5 * np.sqrt(alpha * x / Uinf)

    # Discretized points for y-axis (distance from the plate surface)
    y = np.linspace(0, delta * 2, 100)

    # Velocity Profile: Blasius solution (simplified)
    u = Uinf * (y / delta) ** (1 / 7)

    # Temperature Profile: Assuming linear variation across the thermal boundary layer
    T = Tinf + (Tw - Tinf) * np.exp(-y / delta_t)

    # Plotting the profiles
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # Velocity profile plot
    ax1.plot(u, y, 'b', label="Velocity Profile")
    ax1.set_title('Velocity Profile over Flat Plate')
    ax1.set_xlabel('Velocity (m/s)')
    ax1.set_ylabel('Distance from Plate (m)')
    ax1.grid(True)

    # Temperature profile plot
    ax2.plot(T, y, 'r', label="Temperature Profile")
    ax2.set_title('Temperature Profile over Flat Plate')
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('Distance from Plate (m)')
    ax2.grid(True)

    st.pyplot(fig)
    
# Temperature and Velocity Profile for Flow Inside a Pipe
def pipe_flow_profiles(Vavg, Tw, Tinf, D, Pr):
    nu = 1.5e-5  # Kinematic viscosity (m^2/s) for air at room temp
    alpha = nu / Pr  # Thermal diffusivity

    R = D / 2  # Pipe radius
    r = np.linspace(0, R, 100)  # Radial positions within the pipe

    # Velocity profile (parabolic for laminar flow)
    u = 2 * Vavg * (1 - (r / R) ** 2)

    # Temperature profile (assuming steady-state and laminar flow)
    T = Tw + (Tinf - Tw) * (1 - (r / R) ** 2)

    # Plotting the profiles
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # Velocity profile plot
    ax1.plot(u, r, 'b', label="Velocity Profile")
    ax1.set_title('Velocity Profile in Pipe')
    ax1.set_xlabel('Velocity (m/s)')
    ax1.set_ylabel('Radial Position (m)')
    ax1.grid(True)

    # Temperature profile plot
    ax2.plot(T, r, 'r', label="Temperature Profile")
    ax2.set_title('Temperature Profile in Pipe')
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('Radial Position (m)')
    ax2.grid(True)

    st.pyplot(fig)

# Constants for the methods
def heat_capacity_rate(m_dot, cp):
    return m_dot * cp

def calc_effectiveness(ntu, c_r, flow_type):
    if flow_type == "Parallel Flow":
        return (1 - np.exp(-ntu * (1 + c_r))) / (1 + c_r)
    elif flow_type == "Counter Flow":
        if c_r == 1:
            return ntu / (1 + ntu)
        else:
            return (1 - np.exp(-ntu * (1 - c_r))) / (1 - c_r * np.exp(-ntu * (1 - c_r)))

# ε-NTU Method
def epsilon_ntu_method(mh, mc, cph, cpc, Th_in, Tc_in, U, A, flow_type):
    C_h = heat_capacity_rate(mh, cph)
    C_c = heat_capacity_rate(mc, cpc)
    
    C_min = min(C_h, C_c)
    C_max = max(C_h, C_c)
    
    NTU = U * A / C_min
    C_r = C_min / C_max
    
    epsilon = calc_effectiveness(NTU, C_r, flow_type)
    
    q_max = C_min * (Th_in - Tc_in)
    q = epsilon * q_max
    
    Th_out = Th_in - q / C_h
    Tc_out = Tc_in + q / C_c
    
    return q, Th_out, Tc_out

# LMTD Method
def lmtd_method(Th_in, Tc_in, Th_out, Tc_out, U, A, flow_type):
    if flow_type == "Parallel Flow":
        delta_T1 = Th_in - Tc_in
        delta_T2 = Th_out - Tc_out
    elif flow_type == "Counter Flow":
        delta_T1 = Th_in - Tc_out
        delta_T2 = Th_out - Tc_in

    delta_T_lm = (delta_T1 - delta_T2) / np.log(delta_T1 / delta_T2)
    
    q = U * A * delta_T_lm
    return q

# Plot temperature profiles
def plot_temperature_profiles(Th_in, Tc_in, Th_out_ntu, Tc_out_ntu, Th_out_lmtd, Tc_out_lmtd, flow_type):
    x = np.linspace(0, 1, 100)
    
    if flow_type == "Parallel Flow":
        Th_profile_ntu = Th_in - (Th_in - Th_out_ntu) * x
        Tc_profile_ntu = Tc_in + (Tc_out_ntu - Tc_in) * x

        Th_profile_lmtd = Th_in - (Th_in - Th_out_lmtd) * x
        Tc_profile_lmtd = Tc_in + (Tc_out_lmtd - Tc_in) * x

    elif flow_type == "Counter Flow":
        Th_profile_ntu = Th_in - (Th_in - Th_out_ntu) * x
        Tc_profile_ntu = Tc_out_ntu + (Tc_in - Tc_out_ntu) * x

        Th_profile_lmtd = Th_in - (Th_in - Th_out_lmtd) * x
        Tc_profile_lmtd = Tc_out_lmtd + (Tc_in - Tc_out_lmtd) * x

    fig, ax = plt.subplots(figsize=(10, 5))

    # Plot NTU profiles
    ax.plot(x, Th_profile_ntu, 'r-', label="Hot Fluid (NTU)")
    ax.plot(x, Tc_profile_ntu, 'b-', label="Cold Fluid (NTU)")
    
    # Plot LMTD profiles
    ax.plot(x, Th_profile_lmtd, 'r--', label="Hot Fluid (LMTD)")
    ax.plot(x, Tc_profile_lmtd, 'b--', label="Cold Fluid (LMTD)")
    
    ax.set_title(f"Temperature Profile - {flow_type} (NTU vs. LMTD)")
    ax.set_xlabel("Heat Exchanger Length")
    ax.set_ylabel("Temperature (K)")
    ax.legend()
    ax.grid(True)
    
    st.pyplot(fig)
    
# Streamlit App Interface
st.title('Engineering Analysis Tool')

analysis_type = st.sidebar.selectbox("Choose Analysis", ["Fin Analysis", "Composite Wall Temperature Profile", "Fourier's Law of Heat Conduction", "Transient Cooling with Material Selection", "Planck's and Wien's Law", "3-Body Radiation Problem", "Flow Over Flat Plate", "Flow Inside Pipe", "Heat Exchanger Analysis"])

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
        
elif analysis_type == "Transient Cooling with Material Selection":
    Tinit = st.sidebar.number_input('Enter Initial Temperature (Kelvins)', value=500.0)
    Tf = st.sidebar.number_input('Enter Final Temperature (Kelvins)', value=300.0)
    Tinf = st.sidebar.number_input('Enter Ambient Temperature (Kelvins)', value=293.0)
    h = st.sidebar.number_input('Enter Convective Heat Transfer Coefficient (W/m².K)', value=10.0)
    material = st.sidebar.selectbox("Select Material", ['Iron', 'Copper', 'Aluminum', 'Silver', 'Gold', 'Brass', 'Nickel', 'Lead', 'Zinc', 'Titanium', 'Tungsten', 'Platinum'])
    ro = st.sidebar.number_input('Enter Radius (m)', min_value=0.01, value=0.1)
    
    if st.sidebar.button('Run Transient Cooling Analysis'):
        transient_cooling_with_material_selection(Tinit, Tf, Tinf, h, material)

elif analysis_type == "Planck's and Wien's Law":
    T_values = [st.sidebar.number_input(f'Enter Temperature {i+1} (Kelvins)', value=5000.0) for i in range(3)]
    
    if st.sidebar.button('Run Planck\'s and Wien\'s Law Analysis'):
        planck_wien_law_multi(T_values)
        
elif analysis_type == "3-Body Radiation Problem":
    T1 = st.sidebar.number_input('Temperature of Body 1 (Kelvins)', value=1000.0)
    T2 = st.sidebar.number_input('Temperature of Body 2 (Kelvins)', value=800.0)
    T3 = st.sidebar.number_input('Temperature of Body 3 (Kelvins)', value=600.0)
    e1 = st.sidebar.number_input('Emissivity of Body 1', value=0.9)
    e2 = st.sidebar.number_input('Emissivity of Body 2', value=0.7)
    e3 = st.sidebar.number_input('Emissivity of Body 3', value=0.5)
    A1 = st.sidebar.number_input('Area of Body 1 (m²)', value=1.0)
    A2 = st.sidebar.number_input('Area of Body 2 (m²)', value=1.0)
    A3 = st.sidebar.number_input('Area of Body 3 (m²)', value=1.0)

    if st.sidebar.button('Run 3-Body Radiation Analysis'):
        three_body_radiation(T1, T2, T3, e1, e2, e3, A1, A2, A3)
        
elif analysis_type == "Flow Over Flat Plate":
    Uinf = st.sidebar.number_input('Free-stream Velocity (m/s)', value=10.0)
    Tw = st.sidebar.number_input('Wall Temperature (K)', value=300.0)
    Tinf = st.sidebar.number_input('Ambient Temperature (K)', value=293.0)
    x = st.sidebar.number_input('Distance from Leading Edge (m)', value=0.5)
    Pr = st.sidebar.number_input('Prandtl Number', value=0.71)  # Air has Pr ≈ 0.71

    if st.sidebar.button('Run Flat Plate Flow Analysis'):
        flat_plate_profiles(Uinf, Tw, Tinf, x, Pr)
        
elif analysis_type == "Flow Over Flat Plate":
    Uinf = st.sidebar.number_input('Free-stream Velocity (m/s)', value=10.0)
    Tw = st.sidebar.number_input('Wall Temperature (K)', value=300.0)
    Tinf = st.sidebar.number_input('Ambient Temperature (K)', value=293.0)
    x = st.sidebar.number_input('Distance from Leading Edge (m)', value=0.5)
    Pr = st.sidebar.number_input('Prandtl Number', value=0.71)  # Air has Pr ≈ 0.71

    if st.sidebar.button('Run Flat Plate Flow Analysis'):
        flat_plate_profiles(Uinf, Tw, Tinf, x, Pr)

elif analysis_type == "Flow Inside Pipe":
    Vavg = st.sidebar.number_input('Average Velocity (m/s)', value=2.0)
    Tw = st.sidebar.number_input('Wall Temperature (K)', value=300.0)
    Tinf = st.sidebar.number_input('Fluid Temperature (K)', value=293.0)
    D = st.sidebar.number_input('Pipe Diameter (m)', value=0.1)
    Pr = st.sidebar.number_input('Prandtl Number', value=0.71)  # Air has Pr ≈ 0.71

    if st.sidebar.button('Run Pipe Flow Analysis'):
        pipe_flow_profiles(Vavg, Tw, Tinf, D, Pr)
        
elif analysis_type == "Heat Exchanger Analysis":
    # Select Heat Exchanger Type
    flow_type = st.sidebar.selectbox("Select Flow Type", ["Parallel Flow", "Counter Flow"])
    # Input Parameters
    mh = st.sidebar.number_input('Mass flow rate of hot fluid (kg/s)', value=1.0)
    mc = st.sidebar.number_input('Mass flow rate of cold fluid (kg/s)', value=1.0)
    cph = st.sidebar.number_input('Specific heat capacity of hot fluid (J/kg·K)', value=4186.0)  # Water
    cpc = st.sidebar.number_input('Specific heat capacity of cold fluid (J/kg·K)', value=4186.0)  # Water
    Th_in = st.sidebar.number_input('Inlet temperature of hot fluid (K)', value=370.0)
    Tc_in = st.sidebar.number_input('Inlet temperature of cold fluid (K)', value=300.0)
    U = st.sidebar.number_input('Overall heat transfer coefficient (W/m²·K)', value=500.0)
    A = st.sidebar.number_input('Heat exchanger surface area (m²)', value=10.0)

    # Input for LMTD Method
    Th_out_lmtd = st.sidebar.number_input('Outlet temperature of hot fluid (LMTD) (K)', value=330.0)
    Tc_out_lmtd = st.sidebar.number_input('Outlet temperature of cold fluid (LMTD) (K)', value=320.0)

    if st.sidebar.button("Compare Methods"):
        # ε-NTU Method Calculation
        q_ntu, Th_out_ntu, Tc_out_ntu = epsilon_ntu_method(mh, mc, cph, cpc, Th_in, Tc_in, U, A, flow_type)
        
        # LMTD Method Calculation
        q_lmtd = lmtd_method(Th_in, Tc_in, Th_out_lmtd, Tc_out_lmtd, U, A, flow_type)

        # Display Results
        st.write(f"Heat Transfer Rate (ε-NTU Method): {q_ntu:.2f} W")
        st.write(f"Hot Fluid Outlet Temperature (NTU): {Th_out_ntu:.2f} K")
        st.write(f"Cold Fluid Outlet Temperature (NTU): {Tc_out_ntu:.2f} K")
        st.write(f"Heat Transfer Rate (LMTD Method): {q_lmtd:.2f} W")

        # Plot temperature profiles
        plot_temperature_profiles(Th_in, Tc_in, Th_out_ntu, Tc_out_ntu, Th_out_lmtd, Tc_out_lmtd, flow_type)
