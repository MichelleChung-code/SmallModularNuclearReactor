from math import pi

# Reynolds number
Tavg = (250 + 750) / 2  # 500

rho = 4.1015  # kg/m3
L = 11  # m
mu = 3.8280e-5  # Pa s
D = 3 * .39  # m
mass_flow_rate = 145  # kg/s
A = pi * (D / 2) ** 2
u = mass_flow_rate * A / rho

Re = (rho * u * D) / mu
e = 0.39  # porosity

print("Reynold's Number", Re)
# https://www.hindawi.com/journals/stni/2014/589895/
f_circ = (320 / (Re / (1 - e))) + ((6 / (Re / (1 - e)) ** 0.1))

gc = 1  # dimensional constant m3⋅kg−1⋅s−2

delta_pf = f_circ * (L / D) * (rho * (u ** 2)) / (2 * gc)  # pressure loss

K_in = 0.5
delta_pe_in = K_in * (rho * (u ** 2)) / (2 * gc)  # pressure loss

K_out = 1
delta_pe_out = K_out * (rho * (u ** 2)) / (2 * gc)  # pressure loss

print('Pressure drop across the reactor core:', round((delta_pe_out + delta_pe_in + delta_pf) / 1e3, 2), 'kPa')
