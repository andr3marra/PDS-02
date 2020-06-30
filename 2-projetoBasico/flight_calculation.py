from fluids.atmosphere import ATMOSPHERE_1976
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.constants as sc

#   Import Aerdynamic Data
data = pd.read_csv('~/Documentos/PDS-02/2-projetoBasico/aerodin.dat', skiprows = 1, sep="  ", header=None, engine='python')
data.columns = ("Mach Number", "Cd Rocket - (Unpowered)", "Body Skin Friction", "Nose Pressure Drag", "Base Drag", "Body Pressure Drag", "Interference Drag", "Fin Skin Friction")
#   Simulation
dt          = 0.01
t_end       = 100
n           = int(t_end/dt) + 1
t           = np.linspace( 0, t_end, n)
mass        = np.zeros([n])
acc         = np.zeros([ 3, n])
vel         = np.zeros([ 3, n])
pos         = np.zeros([ 3, n])
teste       = np.zeros([n])
#   Mass
m_0         = 11000             # Dry Mass [g]
#   Engine
mass_f_o    = 2371              # Propelent Mass [g]
mass_flow   = 237.10            # Propelent Mass Flow [g/s]
t_burn      = mass_f_o / mass_flow
thrust      = 500               # Engine Thrust [N]
# Aerodynamic
d           = 152               # mm
A           = sc.constants.pi * (d/1000/2) ** 2

def getDrag(mach):
    index = int(Mach/data["Mach Number"][1])
    mach_part = (mach - data["Mach Number"][index])/(data["Mach Number"][index+1] - data["Mach Number"][index])
    delta_cd = data["Cd Rocket - (Unpowered)"][index+1] - data["Cd Rocket - (Unpowered)"][index]
    Cd = delta_cd * mach_part + data["Cd Rocket - (Unpowered)"][index]
    return Cd

mass[0]     = (m_0 + mass_f_o) / 1000

for i in range( 1,len(t)):                      # Starts on index 1


    # Engine
    if (mass[i-1] > m_0/1000):
        mass[i] = (mass[i-1] - dt * mass_flow/1000)
        T   = thrust 
    else:
        mass[i] = mass[i-1]
        T = 0
    # Atomosphere
    atmo = ATMOSPHERE_1976(pos[0, i])
    # Forces
    W = - mass[i] * atmo.g                                                 # Weight
    # Aerodynamics
    Mach= vel[2, i-1]/atmo.v_sonic
    Cd  =    getDrag(Mach)
    teste[i] = Cd
    D   = -0.5 * atmo.rho * vel[2, i-1] * abs(vel[2, i-1]) * A * Cd           # Drag
    Fz  = W + D + T                                                       # Weight + Drag + Thrust
    acc[:,i] = [0, 0, Fz]/ mass[i]
    # Euler Integration
    vel[:, i] = vel[:, i-1] + acc[:, i] * dt
    pos[:, i] = pos[:, i-1] + vel[:, i] * dt
    # End Simulation
    if (vel[2, i]) < 0:
        i_end = i
        break
#   Torque
    #inercial
    #peso
    #empuxo
    #aerodinamica

# Plot
fig, ax1 = plt.subplots(2)

color = 'tab:red'
ax1[0].set_xlabel('time (s)')
ax1[0].set_ylabel('Acc', color=color)
ax1[0].plot(t[:i_end], acc[2, :i_end], color=color)
#ax1.legend(['label1', 'label2', 'label3'])

ax1[1].set_xlabel('time (s)')
ax1[1].set_ylabel('Vel', color=color)
ax1[1].plot(t[:i_end], vel[2, :i_end], color=color)
#ax1[1].plot(t[:i_end], teste[:i_end], color=color)

ax1[0].tick_params(axis='y', labelcolor=color)

ax2 = ax1[0].twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('h (m)', color=color)  # we already handled the x-label with ax1
ax2.plot(t[:i_end], pos[2, :i_end], color=color)
ax2.tick_params(axis='y', labelcolor=color)

text = "Dry Mass: {} g\nPropelent Mass: {} g\n Max Height: {:.2f} m\n Max Vel {:.2f} m/s\n Max Acc: {:.2f} m/s2".format(m_0, mass_f_o, np.amax(pos[2, :]), np.amax(vel[2, :]), np.amax(acc[2, :]))
fig.text(.74,.32,text,
        horizontalalignment='center')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()