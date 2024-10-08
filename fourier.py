import numpy as np
import matplotlib.pyplot as plt

# Figure properties
width = 6    # Width in inches
height = 5   # Height in inches
alw = 0.75   # AxesLineWidth
fsz = 11     # Fontsize
lw = 1.5     # LineWidth
msz = 8      # MarkerSize

plt.rcParams['lines.linewidth'] = lw
plt.rcParams['lines.markersize'] = msz
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = fsz
plt.rcParams['axes.linewidth'] = alw
plt.rcParams['font.weight'] = 'bold'

# Constants
A = 3   # AREA OF CROSS SECTION (m^2)
K = 2   # THERMAL CONDUCTIVITY OF THE MATERIAL (W/mK)

# Enter parameters
X1 = float(input('Enter the value of X1: '))
X2 = float(input('Enter the value of X2: '))
T1 = float(input('Enter the Temperature of body at X1:\n'))
T2 = float(input('Enter the Temperature of body at X2:\n'))

# Temperature gradient (K/m)
delTX = (T2 - T1) / (X2 - X1)

# Heat transfer rate by conduction (W)
Q = -(K) * (A) * delTX
print(f'THE HEAT TRANSFER RATE BY CONDUCTION IN WATTS IS:\n{Q}')

# Plotting
fig, ax = plt.subplots(figsize=(width, height))
y = np.arange(0, 10.1, 0.1)

# Plotting X1 and X2 positions
if T2 > T1:
    ax.plot(np.ones_like(y) * X1, y, 'bo', label=f'X1={X1}', linewidth=lw)
    ax.plot(np.ones_like(y) * X2, y, 'ro', label=f'X2={X2}', linewidth=lw)
else:
    ax.plot(np.ones_like(y) * X1, y, 'ro', label=f'X1={X1}', linewidth=lw)
    ax.plot(np.ones_like(y) * X2, y, 'bo', label=f'X2={X2}', linewidth=lw)

ax.set_xlim([min([X1, X2]) - 1, max([X1, X2]) + 1])
ax.set_ylim([0, 10])
ax.set_xlabel('Position')
ax.set_ylabel('Y Axis')
ax.set_title('CONDUCTION HEAT TRANSFER (From Red to Blue)')
ax.legend()

plt.show()
