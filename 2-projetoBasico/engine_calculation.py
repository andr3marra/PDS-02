import numpy as np
import scipy.constants as sc
#   Engine
thrust      = 500       #Newtons
P_chamber   = 30        #Bars
L_star      = 1000      #mm
alpha       = 20        #degrees
beta        = 30        #degrees
#   Propelent
Isp         = 240       #s
P_exaust    = 1.025     #Bar
T_chamber   = 2909      #K
o_f         = 3.962     #
gamma       = 1.1574    #
M           = 26.2877   #g/mol

# Calculation
# Propelents
w   = thrust/Isp/sc.constants.pi
w_f = w/(1+o_f)
w_o = w - w_f
R_gas = sc.constants.R*1000/M
#   Throat
T_throat = T_chamber /(1+(gamma-1)/2)
P_throat = P_chamber*(1+(gamma-1)/2)**(-gamma/(gamma-1))
A_throat = (w/P_chamber/100000)*np.sqrt(R_gas*T_throat/gamma)*1000000
print(A_throat)