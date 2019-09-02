"""
MECH4480 Assignment 1
Authors: Patricia Tatel and Tan Thanh Nhan Phan

ABOUT: This Python script is used to investigate the velocity and pressure on 
the ground position

There are four text files to run:
    height_300mm.txt     : Velocity and pressure data at 300mm from the ground, 7 degrees AoA (Far from ground)
    height_15mm.txt      : Velocity and pressure data at 15mm from the ground, 7 degrees AoA (Close to the ground)
    height_70mm.txt      : Velocity and pressure data at 70mm from the ground, 7 degrees AoA ()
    dataground17deg.txt  : Velocity and pressure data at 70mm from the ground, 17 degrees AoA

Change the name of the file at line 53 to run three text files

"""
import matplotlib.pyplot as plt
from math import *
import numpy as np

x = []
y = []
Psi = []
magU = []
U = []
V = []
P = []
Cp = []

#opening raw data which is a txt file
with open('pflowdata.txt') as file:
    for line in file:
        #imported raw txt as 'values' 
        #appended values as floats for math operations
        #eliminated nan values inherted from raw file
        values = line.split()
        if values == []:
            pass
        else:
            x.append(float(values[0]))
            y.append(float(values[1]))
            Psi.append(float(values[2]))
            magU.append(float(values[3]))
            U.append(float(values[4]))
            V.append(float(values[5]))
            P.append(float(values[6]))
            Cp.append(float(values[7]))

groundP = []
groundU = []
groundV = []
groundx =[]
with open('height_70mm.txt') as file:
    for line in file:
        #imported raw txt as 'values' 
        #appended values as floats for math operations
        #eliminated nan values inherted from raw file
        values = line.split()
        if values == []:
            pass
        else:
            groundx.append(float(values[0]))
            groundU.append(float(values[4]))
            groundV.append(float(values[5]))
            groundP.append(float(values[6]))
        

plt.figure()
groundU = np.array(groundU)
groundV = np.array(groundV)
plt.title("Plot of velocity along the ground \n within the design space")
plt.ylabel("Fluid velocity along ground [m/s]")
plt.xlabel("Horizontal position along the design constraints [m]")
plt.grid()
plt.plot(groundx, np.sqrt(groundU**2 + groundV**2), 'r-.')

plt.figure()
plt.title("Ground pressure within the design space")
plt.plot(groundx, groundP, '--.')
plt.ylabel("Pressure [Pa]")
plt.xlabel("Horizontal position along the design constraints [m]")
plt.grid()
plt.show()
