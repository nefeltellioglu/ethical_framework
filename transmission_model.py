#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:18:49 2024

@author: ntelliogluce
"""


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parameters
R0 = 2
I1_infection_duration = 20
I2_infection_duration = 20
initial_pop_size = 1e6
N1_pop_ratio = 0.7
R1_vacc_ratio = 0.2
R2_vacc_ratio = 0.2
vacc_ratio = 0
time_step = 0.0625
final_time = 200

# Initial conditions
I1_0 = 1
I2_0 = 1
R1_0 = 0
R2_0 = 0
S1_0 = initial_pop_size * N1_pop_ratio - R1_vacc_ratio * initial_pop_size * N1_pop_ratio
S2_0 = initial_pop_size * (1 - N1_pop_ratio) - R2_vacc_ratio * initial_pop_size * (1 - N1_pop_ratio)
R1_vaccinated_0 = initial_pop_size * N1_pop_ratio * R1_vacc_ratio
R2_vaccinated_0 = initial_pop_size * (1 - N1_pop_ratio) * R2_vacc_ratio

# Total initial population
N_0 = I1_0 + I2_0 + R1_0 + R2_0 + S1_0 + S2_0

# ODE system
def deriv(y, t, R0, I1_infection_duration, I2_infection_duration, vacc_ratio):
    I1, I2, R1, R2, S1, S2, R1_vaccinated, R2_vaccinated = y
    
    N = I1 + I2 + R1 + R2 + S1 + S2
    
    B_I1_infects_S1 = R0 * I1_infection_duration
    B_I1_infects_S2 = R0 * I2_infection_duration
    B_I2_infects_S1 = R0 * I1_infection_duration
    B_I2_infects_S2 = R0 * I2_infection_duration

    I1_infects_S1 = S1 * I1 * B_I1_infects_S1 / (N )
    I1_infects_S2 = S2 * I1 * B_I1_infects_S2 / (N )
    I2_infects_S1 = S1 * I2 * B_I2_infects_S1 / (N )
    I2_infects_S2 = S2 * I2 * B_I2_infects_S2 / (N )
    
    I1_recovers = I1 / I1_infection_duration
    I2_recovers = I2 / I2_infection_duration
    
    S1_vacc = S1 * vacc_ratio
    S2_vacc = S2 * vacc_ratio
    
    dI1dt = I1_infects_S1 + I2_infects_S1 - I1_recovers
    dI2dt = I1_infects_S2 + I2_infects_S2 - I2_recovers
    dR1dt = I1_recovers
    dR2dt = I2_recovers
    dS1dt = -I1_infects_S1 - I2_infects_S1 - S1_vacc
    dS2dt = -I1_infects_S2 - I2_infects_S2 - S2_vacc
    dR1_vaccinated_dt = S1_vacc
    dR2_vaccinated_dt = S2_vacc
    
    return [dI1dt, dI2dt, dR1dt, dR2dt, dS1dt, dS2dt, dR1_vaccinated_dt, dR2_vaccinated_dt]

# Initial conditions vector
y0 = [I1_0, I2_0, R1_0, R2_0, S1_0, S2_0, R1_vaccinated_0, R2_vaccinated_0]

# Time vector
t = np.arange(0, final_time, time_step)

# Integrate the ODE system
result = odeint(deriv, y0, t, args=(R0, I1_infection_duration, I2_infection_duration, vacc_ratio))

# Plot results
I1, I2, R1, R2, S1, S2, R1_vaccinated, R2_vaccinated = result.T

plt.figure(figsize=(12, 8))
plt.plot(t, I1, label='I1')
plt.plot(t, I2, label='I2')
plt.plot(t, R1, label='R1')
plt.plot(t, R2, label='R2')
plt.plot(t, S1, label='S1')
plt.plot(t, S2, label='S2')
plt.plot(t, R1_vaccinated, label='R1 vaccinated')
plt.plot(t, R2_vaccinated, label='R2 vaccinated')
plt.xlabel('Time (months)')
plt.ylabel('Population')
plt.legend()
plt.title('Disease Spread and Vaccination Dynamics')
plt.grid()
plt.show()
