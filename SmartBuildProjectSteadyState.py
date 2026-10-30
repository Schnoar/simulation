
#General imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import dm4bem

print("running until here")

To = 10                      # Outdoor temperature
h_out = 25               # Convection out
h_in = 8                 # Convection in
l_l = 15                 # length of the outdoor wall in the living room
h_l = 3                  # height of the outdoor wall in the living room
l_b = 5                  # length of the outdoor wall in the bathroom
h_b = 3                  # height of the outdoor wall in the bathroom
S_in = 10                # Surface of the indoor wall betwwen bathroom and living room
S_out_l = l_l * h_l      # Surface of the outdoor wall in the living room
S_out_b = l_b * h_b      # Surface of the outdoor wall in the bathroom
w_in = 0.10              # width of the indoor wall
w_out = 0.15             # width of the outdoor wall
lambda_in = 1.0          # Concrete dense source: https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
lambda_out = 0.5         # Concrete medium source: https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
ACH_window = 1.0         # Air change rate window source: https://cghiaus.github.io/dm4bem_book/tutorials/02_2_1Toy.html#advection
ACH_door = 0.5           # Air change rate door source: https://cghiaus.github.io/dm4bem_book/tutorials/02_2_1Toy.html#advection
V_a_l = 20                  # volume of the air in the thermal zone, m³.
V_a_b = 10                  # volume of the air in the thermal zone, m³.
c_wi = 18216000             # Heat capacity of the indoor wall
c_wo = 	18216000            # Heat capacity of the outdoor wall
rho_indoor = 140         # Density of a concrete wall source: https://www.engineeringtoolbox.com/density-materials-d_1652.html#gsc.tab=0
rho_outdoor = 140        # Density of a concrete wall source: https://www.engineeringtoolbox.com/density-materials-d_1652.html#gsc.tab=0

#HVAC
k_p = 1e06                  # the proportional gain
T_l_sp = 20              # the setpoint temperature of the living room
T_b_sp = 20              # the setpoint temperature of the bathroom



#Matrixes 


#A_Matrix

A = np.zeros([17, 11])       # n° of branches (flows) X n° of nodes(Temperatures)
A[0, 0] = 1                 # branch 0: -> node 0
A[1, 0], A[1, 1] = -1, 1    # branch 1: node 0 -> node 1
A[2, 1], A[2, 2] = -1, 1    # branch 2: node 1 -> node 2
A[3, 2], A[3, 3] = -1, 1    # branch 3: node 2 -> node 3
A[4, 3], A[4, 4] = -1, 1    # branch 4: node 3 -> node 4
A[5, 4], A[5, 5] = -1, 1    # branch 5: node 4 -> node 5
A[6, 5], A[6, 6] = -1, 1    # branch 6: node 5 -> node 6
A[7, 6], A[7, 7] = -1, 1    # branch 7: node 6 -> node 7
A[8, 7], A[8, 8] = -1, 1    # branch 8: node 7 -> node 8
A[9, 8], A[9, 9] = -1, 1    # branch 9: node 8 -> node 9
A[10, 9], A[10, 10] = -1, 1  # branch 10: node 9 -> node 10
A[11, 10] = 1               # branch 11: -> node 10
A[12, 3] = 1               # branch 12: node 3
A[13, 3], A[13, 7] = -1, 1  # branch 13: node 3 -> node 7
A[14, 7] = 1               # branch 14: -> node 7 
A[15, 3] = 1                # branch 15: node 3
A[16, 7] = 1                # branch 16: node 7

# pd.DataFrame(A, index=q, columns=θ)

#G_Matrix

G = np.zeros([17])  
G[0] = h_out * S_out_l                      # Convection outdoor to indoor (wall living rom)
G[1] = lambda_out / w_out * S_out_l         # Conduction  outdoor wall (livining room)
G[2] = lambda_out / w_out * S_out_l         # Conduction  outdoor wall
G[3] = h_in * S_out_l                       # Convection outdoor wall to indoor air
G[4] = h_in * S_in                          # Convection living room to bathroom 
G[5] = lambda_in / w_in * S_in              # Conduction  indoor wall
G[6] = lambda_in / w_in * S_in              # Conduction  indoor wall
G[7] = h_in * S_in                          # Convection living room to bathroom
G[8] = h_out * S_out_b                      # Convection outdoor to indoor (wall bathroom)
G[9] = lambda_out / w_out * S_out_b         # Conduction  outdoor wall (bathroom)
G[10] = lambda_out / w_out * S_out_b        # Conduction  outdoor wall (bathroom)
G[11] = h_in * S_out_b                      # Convection outdoor wall to indoor air (bathroom)
G[12] = ACH_window / 3600 * V_a_l           # Advection window living room
G[13] = ACH_door / 3600 * (V_a_l + V_a_b)   # Advection door living room to bathroom
G[14] = ACH_window / 3600 * V_a_b           # Advection window bathroom
G[15] = k_p                                 # HVAC living room
G[16] = k_p                                 # HVAC bathroom

G = np.diag(G)  # Convert to diagonal matrix
7
#C_Matrix
C = np.zeros([11])
C[0] = rho_outdoor * c_wo * w_in * S_out_l     # Heat capacity of the outdoor wall (living room)
C[1] = rho_indoor * c_wi * w_in * S_in         # Heat capacity of the indoor wall
C[2] = rho_outdoor * c_wo * w_in * S_out_b     # Heat capacity of the outdoor wall (bathroom)

C = np. diag(C)  # Convert to diagonal matrix

# temperature source vector b

b = np.zeros(17)
b[0]  = To
b[1]  = 0
b[2]  = 0
b[3]  = 0
b[4]  = 0
b[5]  = 0
b[6]  = 0
b[7]  = 0
b[8]  = 0
b[9]  = 0
b[10] = 0
b[11] = To
b[12] = To
b[13] = 0
b[14] = To
b[15] = T_l_sp
b[16] = T_b_sp



# flow source vector f

f = np.zeros(11)

θ_steady_To = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)

np.set_printoptions(precision=3)
print('When To = 1°C, the temperatures in steady-state are:', θ_steady_To, '°C')
print(f'The indoor temperature is: {θ_steady_To[-1]:.3f} °C')




print("test")