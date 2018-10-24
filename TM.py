from plotly.offline import iplot, init_notebook_mode 
import numpy as np
import pandas as pd
from scipy import constants
import cmath
init_notebook_mode(connected=True)
import plotly.graph_objs as go
from PyAstronomy.pyaC import zerocross1d


D1, D2, M1, M2, Y1, Y2, Y3, Y4 = 1, 0, 0, 7, 1, 9, 9, 4

c = 3 * np.power(10, 8)
e0 = constants.epsilon_0
u0 = constants.mu_0
impedance_0 = np.sqrt(u0 / e0)

#medium 1
e1r = 1 #relative permitivity
u1r = 1 #relative permiability
n1 = np.sqrt(e1r * u1r) #refractive index
impedance_1 = np.sqrt((u1r * u0)/(e1r * e0))

#medium 2
e2r = (2 + D1 + M2)
u2r = 1
n2 = np.sqrt(e2r * u2r)
impedance_2 = np.sqrt((u2r * u0)/(e2r * e0))

#medium 3
e3r = (Y3 + Y4)/4
u3r = 1
n3 = np.sqrt(e3r * u3r)
impedance_3 = np.sqrt((u3r * u0)/(e3r * e0))

h = 0.01 * (M2 + 10) #thickness of the film
#free_space parameters
#n = np.linspace(start=1, stop=20, num=20)
#omega = (n) * (c /(2 * n2 * h))
omega = (M2 + 1)* (D2 + 1) * np.power(10, 9)
#omega = 0.6302607978880989 * np.power(10, 9)
k0 = (omega/c) #free space wave number

#wave number of each medium

k1 = n1 * k0
k2 = n2 * k0
k3 = n3 * k0

#theta = (np.pi / 6)
theta = np.linspace(0, 0.5 * np.pi, 200)

#tangential wave number:
kt = k0 * np.sin(theta)

#wave number on z direction
kz1 = np.sqrt(np.square(k1) - np.square(kt))
kz2 = np.sqrt(np.square(k2) - np.square(kt))
kz3 = np.sqrt(np.square(k3) - np.square(kt))

cos_theta_1 = kz1 / k1
cos_theta_2 = kz2 / k2
cos_theta_3 = kz3 / k3

phase_real = np.cos(kz2 * h) #real part of the phase
phase_imag = np.sin(kz2 * h) #imaginary part of the phase

phase_real_k3 = np.cos(kz3 * h)
phase_imag_k3 = np.sin(kz3 * h)

# calculation of reflection coefficient:

# ---------------------------------------------------

numerator_real = (impedance_1 * impedance_2 * cos_theta_1 * cos_theta_2) - (impedance_3 * impedance_2 * cos_theta_3 * cos_theta_2)
numerator_imag = (impedance_3 * impedance_1 * cos_theta_3 * cos_theta_1) - (impedance_2 * impedance_2 * cos_theta_2 * cos_theta_2)

denominator_real = (impedance_1 * impedance_2 * cos_theta_1 * cos_theta_2) + (impedance_3 * impedance_2 * cos_theta_3 * cos_theta_2)
denominator_imag = (impedance_3 * impedance_1 * cos_theta_3 * cos_theta_1) + (impedance_2 * impedance_2 * cos_theta_2 * cos_theta_2)

#complex numerator and denominator

complex_numerator = numerator_real * phase_real - 1j * numerator_imag * phase_imag
complex_denominator = denominator_real * phase_real - 1j * denominator_imag * phase_imag

#Reflection Co-eff:

Rvv = (complex_numerator) / (complex_denominator)

Rvv_amplitude = np.absolute(Rvv)
#Rvv_phase = (np.arctan(Rvv.imag / Rvv.real)) / (np.pi)
Rvv_phase = (np.angle(Rvv)) / (np.pi)

#Calculation of Transmission Coefficient

# ---------------------------------------------------

#constant_term = (2 * impedance_3)/(impedance_2 * impedance_1 * cos_theta_2)

#T_numerator_real = (impedance_1 * cos_theta_1) * (np.square(impedance_2 * cos_theta_2 * phase_real) - (impedance_3 * impedance_1 * cos_theta_3 * cos_theta_1 * np.square(phase_imag)))

#T_numerator_imag = (impedance_2 * cos_theta_2) * (np.square(impedance_2 * cos_theta_2) - impedance_3 * impedance_1 * cos_theta_3 * cos_theta_1) * (phase_real * phase_imag)

#complex numerator and denominator

#T_complex_numerator = T_numerator_real - 1j * T_numerator_imag
#T_complex_denominator = complex_denominator

#Transmission Coefficient

#Tvv = constant_term * ((T_complex_numerator) / (T_complex_denominator)) * (phase_real_k3 + 1j * phase_imag_k3)
constant_term = (impedance_3 / (impedance_2 * impedance_1 * cos_theta_2))

T_complex = ((impedance_2 * cos_theta_2) * (1 + Rvv) * phase_real) - 1j * (cos_theta_1 * impedance_1 * (1 - Rvv) * phase_imag) 

T_phase_term = (np.cos(kz3 * h) + 1j * np.sin(kz3 * h))

Tvv = (constant_term) * T_complex * T_phase_term

Tvv_amplitude = np.absolute(Tvv)
Tvv_phase = (np.angle(Tvv)) / (np.pi)

trace1 = go.Scatter(
    x=theta,
    #x=omega,
    y=Rvv_amplitude,
    name='Ref Coeff Amplitude'
)
trace2 = go.Scatter(
    #x=omega,
    x=theta,
    y=Rvv_phase,
    name='Ref Coeff Phase'
)

trace3 = go.Scatter(
    #x=omega,
    x=theta,
    y=Rvv.real,
    name='Ref Real Part'
)

trace4 = go.Scatter(
    #x=omega,
    x=theta,
    y=Rvv.imag,
    name='Ref Imag Part'
)

trace5 = go.Scatter(
    #x=omega,
    x=theta,
    y=Tvv_amplitude,
    name='Transmission Coeff Amplitude'
)
trace6 = go.Scatter(
    #x=omega,
    x=theta,
    y=Tvv_phase,
    name='Transmission Coeff Phase'
)

trace7 = go.Scatter(
    #x=omega,
    x=theta,
    y=Tvv.real,
    name='Transmission Coeff Real'
)
trace8 = go.Scatter(
    #x=omega,
    x=theta,
    y=Tvv.imag,
    name='Transmission Coeff Imag'
)



data = [trace1, trace2, trace3, trace4, trace5, trace6, trace7, trace8]
fig = go.Figure(data=data)
iplot(fig)
