"""
MECH4480 Assignment 1
Authors Patricia Tatel , Tan Thanh Nhan Phan 

About: This script contains function for dicretizing and mapping of aerofoil

Assumptions:
    - Incompressible, inviscid and irrotational flow
Inputs:
    N_vortice          : Number of vortices and constraints along the wing
    U_0                : Freestream velocity (Unit: m/s)
    rho_0              : Freestream density (Unit: kg/m3)
    pin                : Position of mounting point (Unit: m)
    NACAcode           : 4-digit NACA code
    AoA                : Angle of attack (Unit: radians)

Outputs:
    - NACA4421 Aerofoil Geometry
    - Mirror aerofoil plot
    - Vortex strength versus chord position plot
    - Tangential velocity versus chord position plot
    - Total vorticity (Unit: m^2/s)
    - Lift per unit width (Unit: N/m)
    - Total horizontal force Fx (Unit: N/m)
    - Total vertical force Fy (Unit: N/m)
    - Numerical dependence study plot 
    - Number of vortices that meets the threshold (Unit: Dimensionless)
    - Lift coefficient Cl versus angle of attack
"""
#import modules
import numpy as np
import matplotlib.pyplot as plt
from math import * 
import matplotlib.patches as patches

    
            
#functions
def rotate_points(xval, yval, RotM, ref_point):
  """
  ABOUT: This function uses Euler rotation matrix to rotate a coordinate system
  xval, yval about a reference point. 
  
  INPUT:
      xval      : x-coordinate values
      yval      : y-coordinate values
      RotM      : Rotational matrix based on angle of attack AoA
      ref_point : The reference point where the coordinate is shifted
  OUTPUT:
      xnew: new x-coordinate
      ynew: new y-coordinate
  """
  

  #Step 1: Shift the original coordinate system by the offset of
  #the reference point from origin. I.e. shift origin to the point.
  xval = xval - ref_point[0]
  yval = yval - ref_point[1]

  #Step 2: Apply the rotation about this new origin at reference point
  XY = np.matmul(RotM, np.matrix([[xval,yval]]).T)
  
  #Step 3: Shift the coordinate system back to the origin
  xnew = XY[0] + ref_point[0]
  ynew = XY[1] + ref_point[1]


  return xnew, ynew


def NACA_coordinates(NACAcode, chord, N_vortices, alpha, H0, pin):
  """
  ABOUT: A function which is used to generate vortex panels from NACA 4 digit airfoil.
  Each panel has two vortices and one constraint points located at midpoint

  INPUT:
  Input:
    NACAcode     : 4-digit NACA Airfoil Code
    chord        : Chord length of airfoil
    N_vortices   : Number of vortices
    alpha        : Angle of attack - counter clockwise, in radians
    H0           : Distance from ground to the mounting point
    pin          : Location of the mounting point
    
  OUTPUT:
    xupper       : An array of x-coordinate values of upper surface of the wing
    xlower       : An array of x-coordinate values of lower surface of the wing
    yupper       : An array of y-coordinate values of upper surface of the wing
    ylower       : An array of y-coordinate values of lower surface of the wing
    slopeup      : An array of slope angles of upper surface of the wing
    slopedown    : An array of slope angles of lower surface of the wing
  """
  
  #NACA Digits
  M = float(NACAcode[0])/100.
  P = float(NACAcode[1])/10.
  T = float(NACAcode[2:4])/100.      #last 2 digits of Airfoil  
  
  #Set up grid
  Ncalcs = int((N_vortices+1)) 
  beta = np.linspace(0.01, pi, Ncalcs)
  xc = (1.0 - np.cos(beta))/2.0


  #Initialize vectors of locations of vortices and constraint points
  xupper = np.zeros(Ncalcs)
  yupper = np.zeros(Ncalcs)
  xlower = np.zeros(Ncalcs)
  ylower = np.zeros(Ncalcs)
  slopeup = np.zeros(Ncalcs)
  slopedown = np.zeros(Ncalcs)
  theta = np.zeros(Ncalcs)
  
  #constants for thickness distribution and dyt_dx
  a0 = 0.2969
  a1 = -0.126
  a2 = -0.3516
  a3 = 0.2843
  a4 = -0.098
  dxx = 0.005
  
  #Half thickness distribution
  def y_t(x):
    return (T/0.2)*(a0*x**(0.5) + a1*x + a2*x**2.0 + a3*x**3.0 + a4*x**4.0)
  
  def matrix_rotate(alpha):
    """
    ABOUT: A function which is used to compute rotational matrix with initial angle of attack
    INPUT:
          alpha       : Angle of attack (Unit: radians)
    OUTPUT
          RotateMatrix: A rotational matrix evaluated at alpha in anti-clockwise direction    
    """
    RotateMatrix = np.array([[np.cos(alpha), np.sin(alpha)], [-np.sin(alpha), np.cos(alpha)]])
    return RotateMatrix

  RotateMatrix = matrix_rotate(alpha)

                        
   
  #looping of all calculations    
  for i in range(Ncalcs):
      if 0 <= xc[i] < P:
          #front of airfoil - before point of maximum camber pos
          #Camber-line
          yc        = M/(P**2.0)*((2.0*P*xc[i])-xc[i]**2.0) 
          #Gradient
          gradient  = 2.0*M/P**2.0 * (P-xc[i])
          #Theta
          theta[i]  = np.arctan(gradient)

          #Upper surfacee
          xupper[i] = (xc[i] - y_t(xc[i])*np.sin(theta[i]))*chord + dxx
          yupper[i] = (yc + y_t(xc[i])*np.cos(theta[i]))*chord + H0
   			
          #Lower surface
          ylower[i] = (yc - y_t(xc[i])*np.cos(theta[i]))*chord + H0
          xlower[i] = (xc[i] + y_t(xc[i])*np.sin(theta[i]))*chord + dxx
          
          #rotate airfoil
          xupper[i],yupper[i] = rotate_points(xupper[i], yupper[i], RotateMatrix, pin)
          xlower[i],ylower[i] = rotate_points(xlower[i], ylower[i], RotateMatrix, pin)
          
          #derivative and angle
          dyt_dx = (T/0.2)*(0.5*a0*pow(xc[i], -0.5) + a1 + 2*a2*xc[i] + 3*a3*xc[i]**2.0 + 4*a4*xc[i]**3.0) 
          slopeup[i]  = np.arctan(gradient + (dyt_dx*(1.0-theta[i]/2.0) )) - alpha
          slopedown[i] = np.arctan(gradient -(dyt_dx*(1.0-theta[i]/2.0)  )) -alpha
          
      elif P<= xc[i] <= 1.0:
          #back of airfoil - after point of maximum camber pos

          #Camber-line
          yc        = M/(1.0-P)**2.0 *((1-2*P)+2.0*P*xc[i] - xc[i]**2.0)
          #Gradient
          gradient  = 2.0*M/(1-P)**2.0 *(P-xc[i])
          #Theta
          theta[i]  = np.arctan(gradient)
                   
          #Upper surface
          if 0 < i < Ncalcs:
            xupper[i] = (xc[i] - y_t(xc[i])*np.sin(theta[i]))*chord +dxx
            yupper[i] = (yc + y_t(xc[i])*np.cos(theta[i]))*chord + H0

          #Lower surface
          xlower[i] = (xc[i] + y_t(xc[i])*np.sin(theta[i]))*chord + dxx
          ylower[i] = (yc - y_t(xc[i])*np.cos(theta[i]))*chord + H0

          #rotate airfoil
          xupper[i],yupper[i] = rotate_points(xupper[i], yupper[i], RotateMatrix, pin) 
          xlower[i],ylower[i] = rotate_points(xlower[i], ylower[i], RotateMatrix, pin) 
          
          #derivative and angle 
          dyt_dx = (T/0.2)*(0.5*a0*pow(xc[i], -0.5) + a1 + 2*a2*xc[i] + 3*a3*xc[i]**2.0 + 4*a4*xc[i]**3.0) 
          slopeup[i]  = np.arctan(gradient + (dyt_dx*(1.0-theta[i]/2.0)  )) - alpha
          slopedown[i] = np.arctan(gradient - (dyt_dx*(1.0-theta[i]/2.0) )) -alpha

          
          
          if xc[i]==1.0:  
            slopeup[i] = -alpha
            xupper[i] += 0.001
            yupper[i] = (yupper[i-1]+ylower[i-1])/2
            
            
            
 

  return xupper, xlower[1:Ncalcs-1], yupper, ylower[1:Ncalcs-1], slopeup, slopedown[1:Ncalcs-1]


def Const_Vortex(xupper,xlower,yupper, ylower, slopeup, slopedown, AoA):
    """
    ABOUT: A function which is used to compute the vortex and constraint points of the airfoil
    
    INPUT:
        xupper       : An array of x-coordinate values of upper surface of the wing
        xlower       : An array of x-coordinate values of lower surface of the wing
        yupper       : An array of y-coordinate values of upper surface of the wing
        ylower       : An array of y-coordinate values of lower surface of the wing
        slopeup      : An array of slope angles of upper surface of the wing
        slopedown    : An array of slope angles of lower surface of the wing
        AoA          : Angle of attack (Unit: radians)
    OUTPUT:
        xc           : An array of x-coordinate values of constraint points
        yc           : An array of y-coordinate values of constraint points
        xv           : An array of x-coordinate values of vortex points
        yv           : An array of y-coordinate values of vortex points
        slope_angle  : An array of slope angle values of constraint points
    """
    
    #extracting constraints and evaluating constraint angle      
    xc = np.hstack([xupper[1::2], xlower[0::2][::-1]]) 
    yc = np.hstack([yupper[1::2], ylower[0::2][::-1]])
  
    #slope_angle = np.hstack([np.arctan(slopeup[1::2])-AoA, -np.arctan(slopedown[0::2][::-1])-AoA])  
    slope_angle = np.hstack([slopeup[1::2], slopedown[0::2][::-1]])
    #print("angle=", slope_angle)
  
  
    #extracting vortices
    xv = np.hstack([xupper[0::2], xlower[1::2][::-1]])
    yv = np.hstack([yupper[0::2], ylower[1::2][::-1]])
  
    return xc,yc, xv, yv, slope_angle

def Vortex_Strength(x_v ,y_v , x_c, y_c, angle_constraint, AoA,  N_vortices, U, chord):
    """
    ABOUT: A function which is used to compute the vortex strength 
    
    INPUT: 
        x_c               : An array of x-coordinate values of constraint points
        y_c               : An array of y-coordinate values of constraint points
        x_v               : An array of x-coordinate values of vortex points
        y_v               : An array of y-coordinate values of vortex points
        angle_constraint  : An array of slope angle values of constraint points
        AoA               : Angle of attack (Unit: radians)
        N_vortices        : Number of vortices (Unit: Dimensionless)
        U                 : Freestream velocity (Unit: m/s)
        chord             : Chord length of airfoil (Unit: m)
    OUTPUT:
        gamma             : An array of vortex strengths
    """
     
    A = np.zeros((N_vortices, N_vortices)) 
    b = np.zeros(N_vortices)

    #Loop through each constraint points
    for i in range(N_vortices):
        #fluid re-direction
        ac = angle_constraint[i]
			  
        #Create b matrix using the angled uniform flow
        b[i] =  U*np.sin(ac) 
        
        #Loop through each vortices
        for j in range(N_vortices):
            #Compute the square displacement between the ith constraint point and jth vortice
            r2 = pow(x_c[i] - x_v[j], 2.0) + pow(y_c[i] - y_v[j], 2.0)
            r2_mirror = pow(x_c[i] - x_v[j], 2.0) + pow(y_c[i] + y_v[j], 2.0)
           
            #Compute the relative vertical and horizontal velocities 
            Cu  =   1/(2*np.pi*r2) * (y_c[i] - y_v[j]) -  1/(2*np.pi*r2_mirror) * (y_c[i] + y_v[j]) 
            Cv   =  -1/(2*np.pi*r2) * (x_c[i] - x_v[j]) +  1/(2*np.pi*r2_mirror) * (x_c[i] - x_v[j]) 
            
            
            w_ni = Cv*np.cos(ac) - Cu*np.sin(ac);
            #Create A matrix using vertical and horizontal uniform flows
            A[i,j] = w_ni
    #Add the additional linear equation
    gamma = np.dot(np.linalg.inv(A) , b)
    return gamma
    
def Lift_Kutta_Jowkowski(gamma, U_0, rho_0):
    """
    ABOUT: A function which is used to compute the lift based Kutta-Joukowski Theorem
    INPUT: 
        gamma             : An array of vortex strengths
        U_0               : Freestream velocity (Unit: m/s)
        rho_0             : Freestream density (Unit: kg/m3)
    OUTPUT: 
        Lkutta            : Total lift force acting on body based on Kutta-Joukowski Theorem (Unit: N)
    """
    Lkutta = rho_0*U_0*sum(gamma)
    return Lkutta


    
def P_V_distribution(N, U_0, rho_0, alpha, Gamma, x_c, y_c, x_v, y_v):
    """
    ABOUT: A function which is used to compute tangential velocity and pressure distritbution
    INPUT:
        N                 : Number of vortices (Unit: Dimensionless)
        U_0               : Freestream velocity (Unit: m/s)
        rho_0             : Freestream density (Unit: kg/m3)
        alpha             : Angle of attack (Unit: radians)
        gamma             : An array of vortex strengths
        x_c               : An array of x-coordinate values of constraint points
        y_c               : An array of y-coordinate values of constraint points
        x_v               : An array of x-coordinate values of vortex points
        y_v               : An array of y-coordinate values of vortex points
    OUTPUT:
        v_upper           : A list of tangential velocities of upper surface of the wing
        v_lower           : A list of tangential velocities of lower surface of the wing
        P_upper           : A list of pressure of upper surface of the wing
        P_lower           : A list of pressure of lower surface of the wing
    """
    
    #Initialise
    
    P_upper = []
    P_lower = []
    v_upper = []
    v_lower = []
    
    #Velocity components 
    u_0 = U_0*np.cos(alpha)
    v_0 = U_0*np.sin(alpha)
    
    #Trailing edge index
    N_trail = int((N-1)/2)
    
    #Loop through each constraint
    for i in range(N):
        #Re-initialise the velocity components at ith constraint point
        u = 0.0
        v = 0.0
        #Loop through each vortex
        for j in range(N):
          
            r2 = pow(x_c[i] - x_v[j], 2.0) + pow(y_c[i] - y_v[j], 2.0)         
            r2_mirror = pow(x_c[i] - x_v[j], 2.0) + pow(y_c[i] + y_v[j], 2.0)
            
            gammai = Gamma[j]
            
            #Compute the velocity components at ith constraint
            u +=  (y_c[i] - y_v[j]) / (2*pi*r2)*gammai  - (y_c[i] + y_v[j]) / (2*pi*r2_mirror)*gammai
            v +=  -(x_c[i] - x_v[j]) / (2*pi*r2)*gammai + (x_c[i] - x_v[j]) / (2*pi*r2_mirror)*gammai
        v_tangential = np.sqrt(pow(u + u_0, 2) + pow(v + v_0, 2))
        pressure =  -1/2*rho_0*v_tangential**2 
        
        if i <= int((N-1)/2):
            #Upper surface
            v_upper.append(v_tangential)
            P_upper.append(pressure)

        else:
            v_lower.append(v_tangential)
            P_lower.append(pressure)
            
    return v_upper, v_lower, P_upper, P_lower
            
def Lift_per_unit_width(P_upper, P_lower, x_c):
    """
    ABOUT: A function which is used to compute total lift per unit width
    INPUT:
        P_upper           : A list of pressure of upper surface of the wing
        P_lower           : A list of pressure of lower surface of the wing
        x_c               : An array of x-coordinate values of constraint points
    OUTPUT:
        Lift_puw          : Total lift per unit width (Unit: N/m)
    """
    
    #Initalise the values
    pressure_b = 0.0
    pressure_t = 0.0
    #Number of panels on the top and bottom surfave
    N_top = len(P_upper) -1
    N_bot = len(P_lower) -1
    
    #Trapezoidal integration method
    #Calculate the pressure of the top surface of the wing
    for i in range(N_top):
        pressure_t += 0.5*(P_upper[i+1] + P_upper[i])*(x_c[i+1] - x_c[i])
    
    #Calculate the pressure of the bottom surface of the wing
    for i in range(N_bot):
        pressure_b += 0.5*(P_lower[i+1] + P_lower[i])*(x_c[i+1] - x_c[i])
    
    #Total lift per unit width
    Lift_puw = pressure_b - pressure_t#np.trapz(area_b - area_t)
    return Lift_puw

def lift_coefficient(NACAcode,chord, N_vortice, H0, pin, U_0, rho_0):
    """
    ABOUT: A function which is used to compute lift coefficient Cl from a range of 
    angle of attacks 
    
    INPUT:
         NACAcode     : 4-digit NACA Airfoil Code
         chord        : Chord length of airfoil
         N_vortice    : Number of vortices
         H0           : Distance from ground to the mounting point
         pin          : Location of the mounting point
         U_0          : Freestream velocity (Unit: m/s)
         rho_0        : Freestream density (Unit: kg/m3)
    OUTPUT:
        AoA_list      : List of angle of attacks 
        CL_list       : Lift of lift coeffcient
    """
    #Initialise list
    AoA_list =[]
    CL_list = []
    
    AoA_min = np.radians(-15)
    AoA_max = np.radians(15)
    AoA_range = np.linspace(AoA_min, AoA_max, 50, endpoint = True)
    
    for AoA_loop in AoA_range:
        u_0 = U_0 *np.cos(AoA_loop)
        xupper, xlower, yupper, ylower, slopeup, slopedown = NACA_coordinates(NACAcode,chord, N_vortice, AoA_loop, H0, pin)
        xc,yc , xv, yv, anglec = Const_Vortex(xupper,xlower,yupper, ylower, slopeup, slopedown, AoA_loop)
        gamma = Vortex_Strength(xv,yv, xc, yc,anglec, AoA_loop, N_vortice, U_0, chord)
        Lift = Lift_Kutta_Jowkowski(gamma, u_0, rho_0)
        CL = Lift/(0.5*rho_0*u_0**2*chord)
        CL_list.append(CL)
        AoA_list.append(AoA_loop)
        
    return AoA_list, CL_list

def Fy_Fx(P_upper, P_lower, x_v, y_v, angle_constraint, N, alpha):
    """
    ABOUT: A function which is used to compute total vertical and horizontal force 
    acting on a body
    INPUT:
        P_upper           : A list of pressure of upper surface of the wing
        P_lower           : A list of pressure of lower surface of the wing
        x_v               : An array of x-coordinate values of vortex points
        y_v               : An array of y-coordinate values of vortex points
        N                 : Number of vortices
        alpha             : Angle of attack (Unit: radians)
    OUTPUT:
        Fx                : Drag - Force in x direction (Unit: N/m)
        Fy                : Lift - Force in y direction (Unit: N/m)
    """   
    # Initialise the angle of attack rotation matrix (Clockwise rotation).
    RotateMatrix = np.array([[np.cos(alpha), np.sin(alpha)], [-np.sin(alpha), np.cos(alpha)]])
    
    #Join two lists
    P = P_upper + P_lower;
    
    #Initialise the values
    L_total = 0.0 
    D_total = 0.0
    
    # Trailing edge index.
    N_trail = int((N-1)/2);
    
    # Top wing surface.
    for i in range(N_trail+1):
        ds = np.sqrt(pow(x_v[i+1]-x_v[i], 2)+pow(y_v[i+1]-y_v[i], 2));
        L_total  += abs(P[i])*np.cos(angle_constraint[i])*ds;
        
    # Bottom wing surface.
    for i in range(N_trail+1, N-1):
        ds = np.sqrt(pow(x_v[i+1]-x_v[i], 2)+pow(y_v[i+1]-y_v[i], 2));
        L_total  -= abs(P[i])*np.cos(angle_constraint[i])*ds
        
    for i in range(N-1):
        ds = np.sqrt(pow(x_v[i+1]-x_v[i], 2)+pow(y_v[i+1]-y_v[i], 2));
        D_total  += -abs(P[i])*np.sin(angle_constraint[i])*ds;
        
    LiftDrag = np.zeros(2)
    LiftDrag[0] = D_total 
    LiftDrag[1] = L_total
    
    #Perform rotations with angle of attack
    F = np.transpose(np.matmul(RotateMatrix, LiftDrag))
    Fx = F[0]
    Fy = F[1]    
    return Fx, Fy       

def numerical_depedence(NACAcode,chord,AoA, H0, pin, U_0, rho_0):
    """
    ABOUT: A function which is used to compute numerical dependence study to find 
    number of vortices
    INPUT: 
        NACAcode     : 4-digit NACA Airfoil Code
        chord        : Chord length of airfoil
        AoA          : Angle of attack (Unit: radians)
        H0           : Distance from ground to the mounting point
        pin          : Location of the mounting point
        U_0          : Freestream velocity (Unit: m/s)
        rho_0        : Freestream density (Unit: kg/m3)
    OUTPUT:
        None
    """
    #Initialise number of vortices and lift
    N_v = 3
    Lift_p = 0.0
    Lift_c = 0.0
    
    #Data
    Lift = []
    Nv  = []
    u_0 = U_0*np.cos(AoA)
    #v_0 = U_0*np.sin(AoA)
    while True:
        xupper, xlower, yupper, ylower, slopeup, slopedown= NACA_coordinates(NACAcode, chord, N_v, AoA, H0, pin)
        xc,yc , xv, yv, anglec = Const_Vortex(xupper,xlower,yupper, ylower, slopeup, slopedown, AoA)
        
        """Vortex strength"""
        gamma = Vortex_Strength(xv,yv, xc, yc, anglec, AoA,  N_v, U_0, chord)

        """Kutta Lift"""
        Lift_c = Lift_Kutta_Jowkowski(gamma, U_0, rho_0)
        
        if abs(Lift_c - Lift_p) <= 5.0e-4:
            print("Number of vortices that meets requirements = ", N_v)
            break;
        else:
            Lift.append(Lift_c)
            Nv.append(N_v)
            N_v = N_v + 2
            Lift_p = Lift_c
    plt.figure(9)
    plt.plot(Nv, Lift, 'b--')
    plt.title("Numerical dependence study with U_0 = {0}, AoA = {1} with desired Nv = {2}".format(U_0, round(AoA*180/pi,0), N_v))
    plt.ylabel("Lift - Kutta Kuokowski Theorem [Newton]")
    plt.xlabel("Number of vortices (N)")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    """
    Main program
    
    Assumptions:
        1. Incompressible, inviscid and irrotational flow
    Inputs:
            N_vortice          : Number of vortices and constraints along the wing
            U_0                : Freestream velocity (Unit: m/s)
            rho_0              : Freestream density (Unit: kg/m3)
            pin                : Position of mounting point (Unit: m)
            NACAcode           : 4-digit NACA code
            AoA                : Angle of attack (Unit: radians)
    """
    #set up parameters - need to add units
    AoA = radians(7.)			
    H0 = 0.07
    chord = 0.08
    pin = [0.25*chord, H0]
    
    N_vortice = 259 #Odd number
    U_0 = 30.0 
    u_0 = U_0*np.cos(AoA)
    rho_0 = 1.225
    NACAcode = '4421'
    
    
    
    "~~~~~~~ Evaluation ~~~~~~~"
    print("~~~~~~~Running : Airfoil Potential Flow simulation~~~~~~~~~")
    
    """Converting NACA Code to coordinates"""
    xupper, xlower, yupper, ylower, slopeup, slopedown = NACA_coordinates(NACAcode,chord, N_vortice, AoA, H0, pin)
    xc,yc , xv, yv, anglec = Const_Vortex(xupper,xlower,yupper, ylower, slopeup, slopedown, AoA)
    
    plt.figure(2)
    plt.plot(xc,yc, 'rx', label = "Constraints", markersize=5)
    plt.plot(xv, yv, 'bo', label = "Vortices", markersize=5)
    circle1=plt.Circle(pin,.003,color='k',fill=False)
    plt.gcf().gca().add_artist(circle1)
    plt.gcf().gca().add_patch(patches.Rectangle((0, H0 - 0.025), 0.1, 0.05,ls = '--', fill=False ))
    plt.plot(pin[0], pin[1], 'ko', markersize=2)
    plt.legend()
    plt.axis('equal')
    plt.xlabel("Horizontal position (m)")
    plt.ylabel("Vertical position (m)"  )
    plt.title("Schematic of wing geometry and design constraints")
    
    print("Vortex strength plot")
    gamma = Vortex_Strength(xv,yv, xc, yc,anglec, AoA, N_vortice, U_0, chord)
    
    Nt = int((N_vortice-1)/2)
    xv_top = xv[0:Nt+1:1]
    xv_bot = xv[Nt+1:N_vortice:1]
    
    gamma_top = gamma[0:Nt+1:1]
    gamma_bot = gamma[Nt+1:N_vortice:1]
    gamma_sum = sum(gamma)
    plt.figure(3)
    plt.plot(xv_top, gamma_top, 'r-', label = "Upper surface", Markersize=5)
    plt.plot(xv_bot, gamma_bot, 'b', label = "Lower surface", Markersize=5)
    plt.ylabel("Vortex strength $\gamma$  $[m^{2}/s] $")
    plt.xlabel("Position along chord of NACA4421 [m]")
    plt.title("Total circulation $\Gamma = 1.627 m^{2}/s$ at $AoA=7 \deg$")
    plt.grid()
    plt.legend()
    
    """Mirror airfoil"""
    xv_mirror = np.hstack([xv, xv])
    xc_mirror = np.hstack([xc, xc])
    yc_mirror = np.hstack([yc, -1*yc])
    yv_mirror = np.hstack([yv, -1*yv])
    pin_mirror = np.hstack([pin, [pin[0], -1*pin[1]]])
    plt.figure(4)
    plt.plot(xc,yc, 'ro', label = "Constraints no mirror")
    plt.plot(xv, yv, 'bo', label = "Vortices no mirror")
    plt.plot(xc_mirror, yc_mirror, 'r*', label = "Constraints")
    plt.plot(xv_mirror, yv_mirror, 'b*', label = "Vortices")
    plt.plot(pin_mirror[0], pin_mirror[1], 'ko', pin_mirror[2], pin_mirror[3], 'ko')
    plt.legend()
    plt.axis('equal')
    plt.show()
    
    """Pressure and Velocity Distribution"""
    v_upper, v_lower, P_upper, P_lower = P_V_distribution(N_vortice, U_0, rho_0, AoA, gamma, xc, yc, xv, yv)
    xc_top = xc[0:Nt+1:1]
    xc_bot = xc[Nt+1:N_vortice:1]
    
    """Velocity distribution plot"""
    plt.figure(5)
    plt.title("velocity along chord")
    plt.plot(xc_top, v_upper, 'r--', label = "Upper surface")
    plt.plot(xc_bot, v_lower, 'b--', label = "Lower surface")
    plt.title("Tangential Velocity on NACA4421 at $AoA$={} $\deg$".format(round(np.degrees(AoA),0)))
    plt.xlabel("Position across NACA4421 chord [m]")
    plt.ylabel("Velocity [m/s]")
    plt.legend()
    plt.grid()
    plt.show()
    
    """Pressure distribution plot"""
    plt.figure(6)
    xc_top = xc[0:Nt+1:1]
    xc_bot = xc[Nt+1:N_vortice:1]
    plt.plot(xc_top, P_upper, 'r--', label = "Upper surface")
    plt.plot(xc_bot, P_lower, 'b--', label = "Lower surface")
    plt.title("Surface pressure on NACA4421 at $AoA={} \deg$".format(round(np.degrees(AoA),0)))
    plt.xlabel("Position across NACA4421 chord [m]")
    plt.ylabel("Pressure [Pa]")
    plt.grid()
    plt.legend()
    plt.show()
    
    
    """Kutta Lift"""
    lift = Lift_Kutta_Jowkowski(gamma, U_0, rho_0)
    print("Kutta Lift= {:.2f} N".format(lift))
    print("Total vorticity = ", round(gamma_sum, 2), 'm^2/s')
    
    """Lift per unit width"""
    Lift_puw = Lift_per_unit_width(P_upper, P_lower, xc)
    print("Lift per unit width = ", round(Lift_puw,2), "N/m")
    
    """Lift and Drag"""
    Fx , Fy = Fy_Fx(P_upper, P_lower, xv, yv, anglec, N_vortice, AoA)
    print("Total horizontal force Fx = ", round(Fx, 2), "N/m")
    print("Total vertical force Fy = ", round(Fy, 2), "N/m")
    
    AoA_list, CL_list = lift_coefficient(NACAcode,chord, N_vortice, H0, pin, U_0, rho_0)
    
    """Numerical dependence study"""
    numerical_depedence(NACAcode,chord,AoA, H0, pin, U_0, rho_0)
    
    
    """Lift Coefficient plot"""
    plt.figure(7)
    plt.plot(np.degrees(AoA_list), CL_list)
    plt.grid()
    plt.xlabel("Angle of attack (Unit: Degrees)")
    plt.ylabel("Lift coefficient $C_l$ per unit span (Unit: 1/m)")
    plt.title("Lift coeffcient at $U_{\inf}$ = 30.0 m/s, $rho_{0}$ = 1.225 kg/m3")
    plt.show() 

    
    print("Optimum AoA")
    AoA = np.radians(-10)
    AoA_list = []
    lift_list = []
    while AoA < np.radians(60):
        u_0 = 30.0*cos(AoA)
        #v_0 = 30.0*np.sin(AoA)
        xupper, xlower, yupper, ylower, slopeup, slopedown= NACA_coordinates(NACAcode,chord, N_vortice, AoA, H0, pin)
        xc,yc , xv, yv, anglec = Const_Vortex(xupper,xlower,yupper, ylower, slopeup, slopedown, AoA)
        gamma = Vortex_Strength(xv,yv, xc, yc,anglec, AoA, N_vortice, U_0, chord)
        Lkutta = Lift_Kutta_Jowkowski(gamma, u_0, rho_0)
        Cl = Lkutta/(0.5*rho_0*U_0**2*chord)
        #print("Sum of gamma = ", sum(gamma))
        lift_list.append(Cl)
        AoA_list.append(AoA*180/pi)
        AoA += np.radians(1)
    plt.figure(8)
    plt.plot(AoA_list, lift_list)
    plt.title("Lift Coefficient vs. Angle of attack ")
    plt.xlabel("Angle of attack (Degrees)")
    plt.ylabel("Lift coefficient $C_L$ ")
    plt.grid()
    plt.show()    
    
else:
    
    "Importing functions into main"
  
  
