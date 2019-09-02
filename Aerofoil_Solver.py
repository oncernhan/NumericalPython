"""
Assignment 01 - MECH4480 - 18/08/19

Authors: Patricia Tatel and Tan Thanh Nhan Phan


ABOUT: this code imports MECH4480_A01_Vortex.py 
and the functions are evaluated in py script.

HOW TO: to run and save the file along pts
of interest then, use the bash code

pflow.py --job=Aerofoil_Solver.py --out-file==dataout.txt

Required input parameters described in MECH4480_A01_Vortex.py
Both python files must be run within the same folder



"""
#import functions from Aerofoil_PtsFcn
from MECH4480_A01_Vortex import *

#Set model parameters
mdata.name = 'Vortex + Uniform Flow adjacent to a wall'
mdata.dimensions = 2
mdata.rho  = 1.225 


#set up parameters - need to add units
AoA = radians(5.0)			
H0 = 0.07
chord = 0.08
pin = [0.25*chord, H0]
N_vortice = 201 #Odd number
U_0 = 30.0 
u_0 = U_0*np.cos(AoA)
rho_0 = 1.225
NACAcode = '4421'

"~~~~~~~ Evaluation ~~~~~~~"
print("~~~~~~~Running : Airfoil Potential Flow simulation~~~~~~~~~")

  
print("Converting NACA Code to coordinates")
xupper, xlower, yupper, ylower, slopeup, slopedown = NACA_coordinates(NACAcode,chord, N_vortice, AoA, H0, pin)

print("Adding constraints")
xc,yc , xv, yv, anglec = Const_Vortex(xupper,xlower,yupper, ylower, slopeup, slopedown, AoA)

print("Generating vortices")
gamma = Vortex_Strength(xv,yv, xc, yc,anglec, AoA, N_vortice, U_0, chord)

print("Gamma = ", gamma)
plt.figure(5)
plt.plot(xc, gamma, 'y-o')
plt.title("Vortex strength vs Aerofoil location")

print("Vortex strength plot")
Nt = int((N_vortice-1)/2)
xv_top = xv[0:Nt+1:1]
xv_bot = xv[Nt+1:N_vortice:1]

gamma_top = gamma[0:Nt+1:1]
gamma_bot = gamma[Nt+1:N_vortice:1]
  
plt.figure(7)
plt.plot(xv_top, gamma_top, 'r-o', label = "Upper surface")
plt.plot(xv_bot, gamma_bot, 'b-o', label = "Lower surface")
plt.legend()
  
  
"""Ground effects"""
print("Ground effects")
xv_mirror = np.hstack([xv, xv])
xc_mirror = np.hstack([xc, xc])
yc_mirror = np.hstack([yc, -1*yc])
yv_mirror = np.hstack([yv, -1*yv])
pin_mirror = np.hstack([pin, [pin[0], -1*pin[1]]])
plt.figure(6)
plt.plot(xc,yc, 'ro', label = "Constraints no mirror")
plt.plot(xv, yv, 'bo', label = "Vortices no mirror")
plt.plot(xc_mirror, yc_mirror, 'r*', label = "Constraints")
plt.plot(xv_mirror, yv_mirror, 'b*', label = "Vortices")
plt.plot(pin_mirror[0], pin_mirror[1], 'ko', pin_mirror[2], pin_mirror[3], 'ko')
plt.legend()
plt.axis('equal')
  

"""Pressure and Velocity Distribution"""
print("Pressure and Velocity Distribution")
v_upper, v_lower, P_upper, P_lower = P_V_distribution(N_vortice, U_0, rho_0, AoA, gamma, xc, yc, xv, yv)
xc_top = xc[0:Nt+1:1]
xc_bot = xc[Nt+1:N_vortice:1]
  
  
plt.figure(8)
plt.title("velocity along chord")
plt.plot(xc_top, v_upper, 'r-o', label = "Upper surface")
plt.plot(xc_bot, v_lower, 'b-o', label = "Lower surface")
plt.legend()
  
plt.figure(10)
xc_top = xc[0:Nt+1:1]
xc_bot = xc[Nt+1:N_vortice:1]
plt.plot(xc_top, P_upper, 'r-o', label = "Upper surface")
plt.plot(xc_bot, P_lower, 'b-o', label = "Lower surface")

plt.title("Pressure across chord")
plt.legend()
  
plt.show()
  
lift = Lift_Kutta_Jowkowski(gamma, U_0, rho_0)
print("lift=", lift)
#Starage for potential solns
A1 = UniformFlow(U_0, 0.0)

X = []
X.append(A1)

#Evaluate and store vortices
for i in range(N_vortice):
    X.append(Vortex(xv[i],yv[i],Gamma = gamma[i]))


#Evaluate and store vortices
for i in range(N_vortice):
    X.append(Vortex(xv[i],-yv[i],Gamma = -gamma[i]))

# Define how the solution will be visualised. 

visual.xmin= -0.05
visual.xmax = 0.125
visual.ymax = 0.08
visual.ymin = -0.08
visual.subplot = 0
visual.Nx = 300
visual.Ny = 300
plot.psi(levels=50) 
plot.magU(min=0, max=30, levels=30)
plot.P(min=-800, max=0)

#output data 
screen.variables(['Psi','magU','U','V','P','Cp']) 
#along the ground
screen.Lineval([xc[0],xc[-1]], [-0.15,0.15], N=100) 

#if youd like to record along the constraints
"""
for i in range(len(xc)):
    screen.Lineval([xc[i],yc[i]], [xc[i+1],yc[i+1]], N=1)
"""

loc_top = []
loc_bottom = []
for i in range(len(xc)):
    if i == 0:
        loc_top.append([xc[i],yc[i]])        
        loc_bottom.append([xc[i],yc[i]])        
    elif i % 2 == 0:
        loc_bottom.append([xc[i],yc[i]])
    else:
        loc_top.append([xc[i],yc[i]])    






