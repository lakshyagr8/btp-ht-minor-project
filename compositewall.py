# import numpy as np
# import matplotlib.pyplot as plt

# def composite_wall_temperature_profile(Tin, Tout, L1, K1, L2, K2, L3, K3, A):
#     Q = (Tin - Tout) / ((L1 / (K1 * A)) + (L2 / (K2 * A)) + (L3 / (K3 * A)))
#     T1 = Tin - (Q * L1 / (K1 * A))
#     T2 = T1 - (Q * L2 / (K2 * A))
#     T3 = T2 - (Q * L3 / (K3 * A))

#     print(f'TOTAL HEAT TRANSFER THROUGH THE COMPOSITE WALL IN WATTS IS: {Q * 1000}')

#     xmat = np.linspace(0, L1 + L2 + L3, 100)
#     ymat = np.piecewise(
#         xmat,
#         [xmat <= L1, (xmat > L1) & (xmat <= L1 + L2), xmat > L1 + L2],
#         [lambda x: Tin - ((Tin - T1) / L1) * x,
#          lambda x: T1 - ((T1 - T2) / L2) * (x - L1),
#          lambda x: T2 - ((T2 - T3) / L3) * (x - L1 - L2)]
#     )
    
#     plt.figure(figsize=(6, 5))
#     plt.plot(xmat, ymat, 'k')
#     plt.xlabel('x (mm)')
#     plt.ylabel('Temperature (deg C)')
#     plt.grid(True)
#     plt.show()
