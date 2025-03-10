'''
Auxiliary Code for turbulent inflow generation

Generates:
    - Mean inflow profile (velocity, density, temperature)
    - Velocity fluctuation intensity profiles
    - RNG seed values
'''

# from math import exp, log, log10, sin, cos, tan, asin, acos, atan, pi, sqrt, sinh, tanh, cosh
import math # don't overwrite SymPy versions of these math functions
from scipy import interpolate
import sympy
from sympy.parsing.sympy_parser import parse_expr
import numpy as np
from opensbli.core.boundary_conditions.turbulent_inflow.seed_gen import get_seeds
import os
import h5py
import matplotlib.pyplot as plt

#######################################################################################################################
#
# Input initialisation
#
#######################################################################################################################

class input_data(object):
    '''
    Object used to store all the input variables. Initialised with the 'values' list from the set-up file
    '''
    def __init__(self, constants, values):
        input_dictionary = dict(zip(constants, values))
        # varibles from values list
        self.gamma           = float(input_dictionary['gama'])
        self.M               = float(input_dictionary['Minf'])
        self.Re              = float(input_dictionary['Re'])
        self.Pr              = float(input_dictionary['Pr'])
        #self.Twall           = 0.5*(self.Pr)**(1.0/3.0)*(self.M**2.0)*(self.gamma - 1.0) + 1.0
        self.Twall = float(input_dictionary['Twall'])
        self.dt              = float(input_dictionary['dt'])
        self.Nt              = int(input_dictionary['niter'])
        self.Nx              = int(input_dictionary['block0np0'])
        self.Ny              = int(input_dictionary['block0np1'])
        self.Nz              = int(input_dictionary['block0np2'])
        self.SuthT           = float(input_dictionary['SuthT']) # sutherland temperature
        self.RefT            = float(input_dictionary['RefT']) # sutherland reference temperature
        self.betaY           = float(input_dictionary['by']) # grid stretching parameters
        self.betaZ           = float(input_dictionary['bz'])

        self.configuration   = input_dictionary['configuration'] # Type of sidewall: bottom, corner, channel, full, parallel
        self.Re_schlatter    = int(input_dictionary['referenceRe']) # Reynolds number of schlatter target profile
        
        self.Sx              = float(input_dictionary['Lx0']) # Non-dimensionalised grid lengths (x, y, z)
        self.Sy              = float(input_dictionary['Lx1'])
        self.Sz              = float(input_dictionary['Lx2'])
        
        # grid spacing stored in symbolic notation so must be computed to get the float value
        npx, npy, npz        = sympy.symbols(['block0np0','block0np1','block0np2'])
        lx, ly, lz           = sympy.symbols(['Lx0','Lx1','Lx2'])
        exprx                = parse_expr(input_dictionary['Delta0block0'])
        expry                = parse_expr(input_dictionary['Delta1block0'])
        exprz                = parse_expr(input_dictionary['Delta2block0'])
        
        self.dx              = float(exprx.subs([(npx, self.Nx), (lx, self.Sx)])) # Grid spacing
        self.dy              = float(expry.subs([(npy, self.Ny), (ly, self.Sy)]))
        self.dz              = float(exprz.subs([(npz, self.Nz), (lz, self.Sz)]))

        self.wall_normal_v   = False
        
        # Grid arrays are built using the yStretch function
        yGrid = []
        zGrid = []
        deltaT = [self.dt]*self.Nt
        deltaX = [self.dx]*self.Nx
        
        for j in range(self.Ny):
            if self.configuration in ["full", "bottom shock", "channel shock", "parallel"]:
                yGrid.append(zStretch(j, self.betaY, self.Sy, self.Ny))
            else:
                yGrid.append(yStretch(j, self.betaY, self.Sy, self.Ny))

        for k in range(self.Nz):
            if self.configuration=="corner":
                zGrid.append(yStretch(k, self.betaZ, self.Sz, self.Nz))
            else:
                zGrid.append(zStretch(k, self.betaZ, self.Sz, self.Nz))
        
        if self.betaY == 0:
            deltaY = [self.dy]*self.Ny
        else:
            deltaY = []
            for j in range(self.Ny - 1):
                deltaY.append(yGrid[j + 1] - yGrid[j])
            deltaY.append(deltaY[self.Ny - 2])
            
        if self.betaZ == 0:
            deltaZ = [self.dz]*self.Nz
        else:
            deltaZ = []
            for k in range(self.Nz - 1):
                deltaZ.append(zGrid[k + 1] - zGrid[k])
            deltaZ.append(deltaZ[self.Nz - 2])
        
        self.yGrid = yGrid
        self.zGrid = zGrid
        self.deltaT = deltaT
        self.deltaX = deltaX
        self.deltaY = deltaY
        self.deltaZ = deltaZ
        
        self.kapa            = 0.41 # von Karman constant
        self.loglawB         = 5.17 # log-law additive constant
        self.omega           = 0.67 # dynamic viscosity power law coef.

        # Set the viscosity scaling factor
        if "viscosity_scale" in input_dictionary:
            self.vscale = float(input_dictionary['viscosity_scale'])
        else:
            self.vscale = 1.0

        self.NyMean          = 201 # number of y grid points used to build the mean profile
        self.Csuth = self.SuthT/self.RefT

        # Integral length scales
        Lxu                  = 10.0
        Lxv                  = 4.0
        Lxw                  = 4.0
        Lyu                  = 1.5
        Lyv                  = 1.75
        Lyw                  = 1.0
        Lzu                  = 1.5
        Lzv                  = 1.0
        Lzw                  = 1.75

        # Consolidate integral length scales into a single array
        Lscale               = np.array([[Lxu, Lxv, Lxw], [Lyu, Lyv, Lyw], [Lzu, Lzv, Lzw]])
        self.Lscale          = Lscale

    def set_wall_normal_v(self, condition):
        self.wall_normal_v = condition
        return
            
def yStretch(j, beta, Sy, Ny):
    '''
    asymmetric stretching function
    '''
    if float(beta) == 0.0:
        # Exception for when stetch factor is zero (no stetching)
        return Sy*j/(Ny - 1)
        
    else:
        # Implements the hyperbolic sine function used in openSBLI

        return Sy * math.sinh(beta*j /(Ny - 1))/math.sinh(beta)
        # return Sy*0.5*(1 - math.tanh(beta*(1.0 - 2.0*j/(2*Ny - 1.0)))/math.tanh(beta)) 

def zStretch(k, beta, Sz, Nz):
    '''
    symmetric stretching function
    '''
    if float(beta) == 0.0:
        # Exception for when stetch factor is zero (no stetching)
        return Sz*k/(Nz - 1)
        
    else:
        # Implements the hyperbolic sine function used in openSBLI

        return Sz*0.5*(1 - math.tanh(beta*(1.0 - 2.0*k/(Nz - 1.0)))/math.tanh(beta))    
        
#######################################################################################################################
#
# Auxiliary functions
#
#######################################################################################################################


def friction_velocity(inputs, Umean, Tmean, Dmean, grid='y'):
    '''
    Function which calculates the friction (shear) velocity profile of the mean profile - used to computed wall-normal values
    '''
    Ny = inputs.Ny
    yGrid = inputs.yGrid
    zGrid = inputs.zGrid
    Csuth = inputs.Csuth
    vscale = inputs.vscale
    
    u0 = 0
    u1 = Umean[1]
    u2 = Umean[2]
    u3 = Umean[3]
    u4 = Umean[4]

    if grid=='y':
        dy0 = yGrid[1] - yGrid[0]
    else:
        dy0 = zGrid[1] - zGrid[0]

    # compute wall velocity gradient using the 3-order forward difference approximation of the first derivative
    # dudy = (-11.0*u0/6.0 + 3.0*u1 - 1.5*u2 + u3/3.0)/dy0
    dudy = (-25.0*u0/12.0 + 4.0*u1 - 3.0*u2 + 4.0*u3/3.0 - u4/4.0)/dy0

    mu = vscale*(Tmean[0]**1.5)*(1 + Csuth)/(Tmean[0] + Csuth)
    # wall kinematic viscosity
    kinVisc = mu/Dmean[0]
    uTau = (kinVisc*dudy)**0.5
     
    return uTau, kinVisc

def plus_values(y, u, uTau, kinVisc, Ny):
    '''
    Function which calculates the wall coordinate system and the dimenstionless velocity profile
    
    The input u is any velocity profile (can be mean, fluctuating or total)
    '''
    # predefine arrays
    yPlus = np.zeros(Ny)
    uPlus = np.zeros(Ny)
    
    # Compute values
    for j in range(Ny):
        yPlus[j] = y[j]*uTau/kinVisc
        uPlus[j] = u[j]/uTau
        
    return yPlus, uPlus

def plus_invert(y_plus, u_plus, uTau, kinVisc, Ny):
    '''
    Function which converts wall-normalised values into physical values
    '''
    # predefine arrays
    y = np.zeros(Ny)
    u = np.zeros(Ny)
    
    # Compute values
    for j in range(Ny):
        y[j] = y_plus[j]*kinVisc/uTau
        u[j] = u_plus[j]*uTau
        
    return y, u
       
#######################################################################################################################
#
# Mean profile
#
#######################################################################################################################

def VD_profile(inputs, u_VD_eplus, option=[False]):
    '''
    Function that computes the velocity profile in van Driest coordinates
    '''
    if option[0]:
        Re = option[1]
    else:
        Re = inputs.Re
    NyMean = inputs.NyMean
    kapa = inputs.kapa
    loglawB = inputs.loglawB
    
    # Initialise arrays
    u_VD_plus = [0.0]*NyMean # wall normalised streamwise velocity in VD space
    u_VD = [0.0]*NyMean # velocity in VD space
    xi_plus = [0.0]*NyMean
    y_plus = [0.0]*NyMean # wall normal coordinates in VD space
    y_VD = [0.0]*NyMean # coordinates in VD space
    
    # Lay out VD velocity distribution and compute xi values 
    for i in range(NyMean):
        u_VD_plus[i] = u_VD_eplus*i/(NyMean - 1)
        tmp = u_VD_plus[i]
        
        xi_plus[i] = tmp + math.exp(-kapa*loglawB)*(math.exp(kapa*tmp) - 1 - kapa*tmp - 0.5*(kapa*tmp)**2 - ((kapa*tmp)**3)/6)
    
    # Compute the new coordinate system 
    xi_eplus = xi_plus[NyMean - 1]
    
    # Compute the new coordinate system (wall normalised)
    for i in range(NyMean - 1):
        y_plus[i] = xi_plus[i]*xi_eplus/(xi_eplus - xi_plus[i]) 
    
    y_plus[NyMean - 1] = 1000000 # final coordinate is a sufficiently large value representing far-field
    
    
    # coordinates in VD space - now not wall normalised    
    for i in range(NyMean):
        y_VD[i] = y_plus[i]*u_VD_eplus/Re
    
    # Normalised boundary layer thickenss
    alpha = math.exp(2*Re/(690 + 1.5*Re)) - 1
    delta0 = u_VD_eplus*xi_eplus/(alpha*Re)
    
    # Velocity profile in VD space (normalised by free stream velocity)
    u_VD[0] = 0.0
    for i in range(1, NyMean):
        eta = y_VD[i]/delta0
        try:
            f = math.exp(-3*(math.exp(eta**(1.0/kapa)) - 1 )) 
            u_VD[i] = 1 - f + (u_VD_plus[i]/u_VD_eplus)*f
        except:
            u_VD[i] = 1.0
        
        
    return u_VD, y_VD          
                
def physical_profile(inputs, u_VD_eplus, u_VD, y_VD):
    '''
    Function that computes the steamwise velocity, temperature and density profiles based on VD velocity profile
    '''
    NyMean = inputs.NyMean
    gamma = inputs.gamma
    M = inputs.M
    Twall = inputs.Twall
    Csuth = inputs.Csuth
    vscale = inputs.vscale
    
    # Predefine variables with empty arrays
    u = [0]*NyMean
    T = [0]*NyMean
    rho = [0]*NyMean
    y = [0]*NyMean
    
    # Method for generating the mean temperature profile
    Taw = 1 + 0.5*(gamma - 1)*(M**2)
    Tw = Twall
    a1 = Tw
    b1 = 1 + 0.5*(gamma - 1)*(M**2) - Tw
    c1 = -0.5*(gamma - 1)*(M**2) 
    
    # Find the edge velociy in + units
    K = math.sqrt(-Tw/c1)*(math.asin(b1/math.sqrt(b1**2 - 4*a1*c1)) - math.asin((2*c1 + b1)/math.sqrt(b1**2 - 4*a1*c1)))
    u_eplus = u_VD_eplus/K 
    
    # Iterate through boundary layer to find velocity, temperature and density
    for i in range(NyMean):
        f = u_VD[i]*u_VD_eplus/u_eplus 
        K = math.asin(b1/math.sqrt(b1**2 - 4*a1*c1)) - f*math.sqrt(-c1/Tw) 
        u[i] = (math.sqrt(b1**2 - 4*a1*c1)*math.sin(K) - b1)/(2*c1)
        T[i] = a1 + b1*u[i] + c1*(u[i]**2) 
        rho[i] = 1/T[i]
    
    # Method to convert grid from VD space into real space - uses either Sutherland law or Power law   
    for i in range(NyMean):
    #   power  law  version
    #   y[i] = u_eplus*(T[0]**(omega + 1))*y_VD[i]/u_VD_eplus
    #   Sutherland  law  version 
        y[i] = u_eplus*(vscale*(T[0]**(5.0/2))*((1 + Csuth)/(T[0] + Csuth)))*y_VD[i]/u_VD_eplus 
    
    return u, T, rho, y


def deltaVD(u_VD, y_VD, NyMean):
    '''
    Function which finds the VD displacement thickness by integrating u_VD
    
    '''
    delta1_VD = 0.0
    for i in range(NyMean - 1):
        delta1_VD = delta1_VD + 0.5*(2 - u_VD[i] - u_VD[i + 1])*(y_VD[i + 1] - y_VD[i])    
    
    return delta1_VD     
        
def mean_inflow_profile(inputs, direction='y'):
    '''
    Function that calculates the mean velocity, temperature & density profiles at the inflow
    
    The method first must find the correct u_VD and y_VD profiles which produce a delataVD value of 1.0
    
    The technique used here has upper and lower bound guesses of u_VD_eplus. A new guess, half-way between, is then analysed to find a value of deltaVD. If this is higher/lower than 1.0 then te guessed value of u_VD_eplus becomes the new upper/lower bound value. This process is repeaded until it becomes sufficiently converged
    '''
    
    Re = inputs.Re
    NyMean = inputs.NyMean
    Ny = inputs.Ny
    Nz = inputs.Nz
    yGrid = inputs.yGrid
    zGrid = inputs.zGrid
    
    # First make sure that the first upper/lower-bound values produce deltaVD of respectively greater/less than than 1.0
    
    # upper and lower bound values for guessing
    upper = 30.0
    lower = 20.0
    del_x = 0.25
    
    while True:
        
        # Compute the profiles and deltaVD value
        u_upper, y_upper = VD_profile(inputs, upper)
        d_upper = deltaVD(u_upper, y_upper, NyMean)
        
        u_lower, y_lower = VD_profile(inputs, lower)
        d_lower = deltaVD(u_lower, y_lower, NyMean)
        
        flag = False
        
        if d_upper < 1.0:
            upper *= 2
            flag = True
            
        if d_lower > 1.0:
            lower *= 0.5
            flag = True
            
        if not flag:
            break
            
    # Margin of error for converged value of deltaVD
    e = 0.0001
    
    # Implementation of iterative method
    while True:
        # Guess value of u_VD_eplus
        uGuess = 0.5*(upper + lower)
        # compute VD profiles
        u_VD, y_VD = VD_profile(inputs, uGuess)
        d = deltaVD(u_VD, y_VD, NyMean)
        
        # diference between required and actual value
        check = d - 1.0
        
        if (abs(check) < e):
            # if within error of margin, process has convereged sufficiently
            break
        
        # Find new upper/lower bound value    
        if (check > 0):
            upper = uGuess
        else:
            lower = uGuess
    
    # From the VD profiles, the physical profiles can be computed   
    u_VD_eplus = uGuess
    u, T, rho, y = physical_profile(inputs, u_VD_eplus, u_VD, y_VD)        

    if inputs.wall_normal_v:
        v = wall_normal_new(inputs, y, u, T, u_VD_eplus)
    else:
        v = np.zeros(NyMean)
    
    '''
    The final step is to use spline interpolation to convert the u, v, T & rho profiles to match the grid coordinates used in the main method.
    '''
    
    # Create splines of each of the profiles using the interpolate.splrep (spline representation) function from the scipy library
    splineU = interpolate.splrep(y, u, s=0)
    splineV = interpolate.splrep(y, v, s=0)
    splineT = interpolate.splrep(y, T, s=0)
    splineD = interpolate.splrep(y, rho, s=0)
    
    # From these spline representations, the required mean profiles can be found using in the interpolate.splev (spline evaluate) function
    if direction=='y':
        Umean = interpolate.splev(yGrid, splineU, der=0)
        Tmean = interpolate.splev(yGrid, splineT, der=0)
        Dmean = interpolate.splev(yGrid, splineD, der=0)
        Vmean = interpolate.splev(yGrid, splineV, der=0)
    else:
        Umean = interpolate.splev(zGrid, splineU, der=0)
        Tmean = interpolate.splev(zGrid, splineT, der=0)
        Dmean = interpolate.splev(zGrid, splineD, der=0)
        Vmean = interpolate.splev(zGrid, splineV, der=0)
    
    return Umean, Vmean, Tmean, Dmean

def upstream_profile(inputs, dx, uvdeplus, u, T, y):

    NyMean = inputs.NyMean
    
    # find incompressible skin friction from slope of velocity profile
    u00, u01, u02, u03 = u[0], u[1], u[2], u[3]
    dy0 = y[1] - y[0]
    dudy = (-11.0*u00/6.0 + 3.0*u01 - 1.5*u02 + u03/3.0)/dy0
    cf1 = 2*dudy/inputs.Re

    cf1 = 0.00333
    #print cf1

    # Use power law from White (Visc. Fluid Flow, 1991) to find Cf at point (2)
    Re1 = (cf1/0.020)**(-6)
    Rex1 = (Re1/0.16)**(7.0/6.0)
    x1 = Rex1/inputs.Re
    Rex2 = Rex1*(1.0 + dx/x1)
    Re2 = 0.16*(Rex2**(6.0/7.0))
    cf2 = 0.02*(Re2**(-1.0/6.0))#*inputs.Re

    # Get U_{vd,e}^+ from ratio of skin friction at point (1) and (2)
    uvdeplus2 = uvdeplus*math.sqrt(cf1/cf2)

    #print x1

    # Get delta_{vd} at point (2) based on 1/5ths power law
    delta_vd2 = (1.0 + dx/x1)*(Rex1/Rex2)**(0.2)
    Re_vd2 = inputs.Re#*delta_vd2

    # print uvdeplus2/uvdeplus
    # print delta_vd2

    #print uvdeplus2
    #print delta_vd2
    
    # From this, the downstream physical profile can be computed
    uvd2, yvd2 = VD_profile(inputs, uvdeplus2, option=[True, Re_vd2])
    u2, T2, rho2, y2 = physical_profile(inputs, uvdeplus2, uvd2, yvd2)
    
    # print y
    # print y2

    ### interpolate on to original grid
    spline1 = interpolate.splrep(y2, u2, s=0)
    spline2 = interpolate.splrep(y2, T2, s=0)
    
    u2 = interpolate.splev(y, spline1, der=0)
    T2 = interpolate.splev(y, spline2, der=0)

    return u2, T2

def wall_normal_new(inputs, y, u, T, uvdeplus):
    # calculate v-velocity profile at the inflow
    # downstream profiles are generated and used to find  du/dx(y) at x = 0
    # the v distribution is found by integrating the continuity equation

    NyMean = inputs.NyMean
    rhou = np.zeros((3, NyMean))
    
    # generate upstream profiles to find rho*u at different x locations
    for i, dx in enumerate([1.0, 2.0, 3.0]):
        u2, T2 = upstream_profile(inputs, dx, uvdeplus, u, T, y)

        rhou[i, :] = u2/T2

    ru0 = np.array(u)/np.array(T)
    ru1, ru2, ru3 = rhou

    dx = 1.0
    v = np.zeros(NyMean)

    # find du/dx at x=0 using a 4th order forward difference
    dudx = ((-11.0/6)*ru0 + (3)*ru1 - (1.5)*ru2 + (1.0/3)*ru3)/dx

    # integration of continuity equation to find v(y)
    for i in range(1, NyMean):
        v[i] = v[i - 1] - 0.5*T[i]*(dudx[i] + dudx[i - 1])*(y[i] - y[i - 1])

    # write to file
    f = h5py.File('dudx.h5', 'w')
    f.create_dataset("y", data=y)
    f.create_dataset("dudx", data=dudx)
    f.create_dataset("v", data=v)
    f.close

    #print v

    return np.array(v)

def wall_normal_profile(inputs, y, u, T, u_VD_eplus):
    '''
    Calculate V(y) distribution
    
    The final distribtion to calculate is the wall-normal velocity profile.
    
    First, a new (downstream) profile must be generated with new values of u_VD_eplus and Re_VD
    '''
    NyMean = inputs.NyMean
    del_x = 0.25
     
    # Implemented here is an empirical method to find the new VD values
    Cf1 = 2.0/(u_VD_eplus**2)
    Re1 = (Cf1/0.020)**(-6)
    ReX1 = (Re1/0.16)**(7.0/6.0)
    X1 = ReX1/inputs.Re
    ReX2 = ReX1*(1.0 + del_x/X1)
    Re2 = 0.16*(ReX2**(6.0/7.0))
    Cf2 = 0.02*(Re2**(-1.0/6.0))

    u_VD_eplus2 = math.sqrt(2.0/Cf2)
    Re_VD2 = ReX2/(X1 + del_x)
    
    # From this, the downstream physical profile can be computed
    u_VD2, y_VD2 = VD_profile(inputs, u_VD_eplus2, option=[True, Re_VD2])
    u2, T2, rho2, y2 = physical_profile(inputs, u_VD_eplus2, u_VD2, y_VD2)
    
    '''
    Now the method uses the continuity equation to find the distribution of wall-normal velocity
    ''' 
   
    # Predefine arrays
    Uint = [0.0]*NyMean
    v = [0.0]*NyMean
    Tint = [0.0]*NyMean
    dudx = [0.0]*NyMean
    Tint[0] = T2[0]
    
    # Interpolation between the two profiles
    for i in range(1, NyMean):
        # find  two  nearest  y2's  of  y_j
        level = y[i]
        n = 0
        for jj in range(NyMean - 1):
            dif1 = level - y2[jj] 
            dif2 = level - y2[jj + 1]
            if (dif1*dif2 < 0.0):
                n = jj 
        # linear  interpolation 
        dif1 = level - y2[n]
        dif2 = y2[n + 1] - level 
        dist = dif1 + dif2
        Uint[i] = u2[n]*(dif2/dist) + u2[n + 1]*(dif1/dist) 
        Tint[i] = T2[n]*(dif2/dist) + T2[n + 1]*(dif1/dist)
        
    # Find the distribution of d(rho*u)/dx
    for i in range(NyMean):
        dudx[i]=(Uint[i]/Tint[i]-u[i]/T[i])/del_x
    
    # Continuity equation is integrated here to find the distribution of v(y)
    for i in range(1, NyMean):
        v[i] = v[i - 1] - 0.5*T[i]*(dudx[i] + dudx[i - 1])*(y[i] - y[i - 1])

    return np.array(v)

def centerline_damping(inputs, Umean_y, Umean_z):
    ## Produces a damping matrix to reduce the wall normal velocities towards the centreline (where there are opposing walls)
    Ny = inputs.Ny
    Nz = inputs.Nz

    ygrid = inputs.yGrid
    zgrid = inputs.zGrid

    ney = int(np.interp(0.99, Umean_y, np.linspace(0,Ny-1,Ny)))
    nez = int(np.interp(0.99, Umean_z, np.linspace(0,Nz-1,Nz)))

    bley = np.interp(0.99, Umean_y, ygrid)
    blez = np.interp(0.99, Umean_z, zgrid)

    dampy = np.ones(Ny)
    dampz = np.ones(Nz)

    for j in range(ney,int(0.5*Ny)):
        dampy[j] = (ygrid[j] - 0.5*max(ygrid))/(bley - 0.5*max(ygrid))

    for k in range(nez,int(0.5*Nz)):
        dampz[k] = (zgrid[k] - 0.5*max(zgrid))/(blez - 0.5*max(zgrid))


    return dampy, dampz


def blend_mean_profile(inputs):

    Re = inputs.Re
    Re_theta = inputs.Re_schlatter
    Ny = inputs.Ny
    Nz = inputs.Nz
    configuration = inputs.configuration

    Umean_y, Vmean_y, Tmean_y, Dmean_y = mean_inflow_profile(inputs)
    Umean_z, Wmean_z, Tmean_z, Dmean_z = mean_inflow_profile(inputs, direction='z')

    dampy, dampz = centerline_damping(inputs, Umean_y, Umean_z)

    if any(conf in configuration for conf in ["full", "parallel"]):
        ### check for y mirroring
        mid_height = int(np.ceil(Ny/2))
        Vmean_y = Vmean_y*dampy
        for j in range(mid_height, Ny):
            Umean_y[j] = Umean_y[Ny - j - 1]
            Vmean_y[j] = -Vmean_y[Ny - j - 1]
            Tmean_y[j] = Tmean_y[Ny - j - 1]
            Dmean_y[j] = Dmean_y[Ny - j - 1]
    if any(conf in configuration for conf in ["bottom", "parallel"]):
        ### cases with no corner blending
        Umean = np.repeat(Umean_y[:,np.newaxis], Nz, axis=1)
        Tmean = np.repeat(Tmean_y[:,np.newaxis], Nz, axis=1)
        Dmean = np.repeat(Dmean_y[:,np.newaxis], Nz, axis=1)
        Vmean = np.repeat(Vmean_y[:,np.newaxis], Nz, axis=1)
        Wmean = np.zeros((Ny,Nz))
    else:
        ### cases with corner blending
        Umean = np.zeros((Ny,Nz))
        Vmean = np.zeros((Ny,Nz))
        Wmean = np.zeros((Ny,Nz))
        Tmean = np.zeros((Ny,Nz))
        Dmean = np.zeros((Ny,Nz))
        if any(conf in configuration for conf in ["full", "channel"]):
            ### cases with z mirroring
            mid_span = int(np.ceil(Nz/2))
            Wmean_z = Wmean_z*dampz
            for k in range(mid_span, Nz):
                Umean_z[k] = Umean_z[Nz - k - 1]
                Wmean_z[k] = -Wmean_z[Nz - k - 1]
                Tmean_z[k] = Tmean_z[Nz - k - 1]
                Dmean_z[k] = Dmean_z[Nz - k - 1]

        for z in range(Nz):
            for y in range(Ny):
                ### corner blending
                Tw = inputs.Twall
                a1 = Tw
                b1 = 1 + 0.5*(0.4)*(inputs.M**2) - Tw
                c1 = -0.5*(0.4)*(inputs.M**2) 
                Umean[y,z] = Umean_z[z]*Umean_y[y]
                Vmean[y,z] = Umean_z[z]*Vmean_y[y]
                Wmean[y,z] = Umean_y[y]*Wmean_z[z]
                Tmean[y,z] = a1 + b1*Umean[y,z] + c1*(Umean[y,z]**2) 
                Dmean[y,z] = 1.0/Tmean[y,z]

    # Vmean = np.zeros((Ny,Nz))
    # Wmean = np.zeros((Ny,Nz))

    return Umean, Vmean, Wmean, Tmean, Dmean


#######################################################################################################################
#
# Fluctuatating & combined profile
#
#######################################################################################################################            
    

def import_stats(Re):
    '''
    Reads target profile data from reference file and returns required variables.
    '''

    names = ['y', 'y_plus', 'U_plus', 'urms_plus', 'vrms_plus', 'wrms_plus', 'uv_plus']

    hf = h5py.File('reference.h5', 'r')
    group = hf[str(Re)]

    ret = []
    for n in names:
        ret.append(np.array(group[n]))

    return ret
    # return y, y_plus, U_plus, urms_plus, vrms_plus, wrms_plus, uv_plus
    
def process_fluctuation_input(inputs, Umean, Tmean, Dmean, direction='y'):
    '''
    Function used to process the input data into physical fluctuations
    '''
    Re_theta = inputs.Re_schlatter
    Re = inputs.Re

    if direction=='y':
        yGrid = inputs.yGrid
        Ny = inputs.Ny
    else:
        yGrid = inputs.zGrid
        Ny = inputs.Nz
    
    # Gather data from schlatted DNS files
    fname = 'schlatter_data_{}.txt'.format(Re_theta)
    [y_schlatter, y_plus_schlatter, U_plus_schlatter, urms_plus_schlatter, vrms_plus_schlatter, wrms_plus_schlatter, uv_plus_schlatter] = import_stats(Re_theta)
    
    # Get wall-normalised data from mean profile
    uTau, kinVisc = friction_velocity(inputs, Umean, Tmean, Dmean, grid=direction)
    yPlus, uPlus = plus_values(yGrid, Umean, uTau, kinVisc, Ny)
    yPlus *= Re**0.5

    # Code block to extrapolate the reference data for large wall-normal distances
    if yPlus[-1] > y_plus_schlatter[-1]:
        delta = 100*(y_plus_schlatter[-1] - y_plus_schlatter[-2])
        u_final = 1E-06

        y_extend = [y_plus_schlatter[-1]]
        n = 0
        while y_extend[-1] < yPlus[-1]:
            y_last = y_extend[-1]
            y_extend.append(y_last + delta)
            n += 1

        y_plus_schlatter = np.append(y_plus_schlatter, y_extend[1:])
        urms_plus_schlatter = np.append(urms_plus_schlatter, np.linspace(urms_plus_schlatter[-1], u_final, n+1)[1:])
        vrms_plus_schlatter = np.append(vrms_plus_schlatter, np.linspace(vrms_plus_schlatter[-1], u_final, n+1)[1:])
        wrms_plus_schlatter = np.append(wrms_plus_schlatter, np.linspace(wrms_plus_schlatter[-1], u_final, n+1)[1:])
        uv_plus_schlatter = np.append(uv_plus_schlatter, np.linspace(uv_plus_schlatter[-1], u_final, n+1)[1:])
    
    # Use splineinterpolation to match data to mean profile coordinates
    spline1 = interpolate.splrep(y_plus_schlatter, urms_plus_schlatter, s=0)
    spline2 = interpolate.splrep(y_plus_schlatter, vrms_plus_schlatter, s=0)
    spline3 = interpolate.splrep(y_plus_schlatter, wrms_plus_schlatter, s=0)
    spline4 = interpolate.splrep(y_plus_schlatter, uv_plus_schlatter, s=0)
    
    urms_plus = interpolate.splev(yPlus, spline1, der=0)
    vrms_plus = interpolate.splev(yPlus, spline2, der=0)
    wrms_plus = interpolate.splev(yPlus, spline3, der=0)
    uv_plus = interpolate.splev(yPlus, spline4, der=0)
    
    # density scaling to convert to compressible
    for i in range(Ny):
        dratio = ((Dmean[0]/Dmean[i])/Re)**0.5
        urms_plus[i] *= dratio
        vrms_plus[i] *= dratio
        wrms_plus[i] *= dratio
        uv_plus[i] *= dratio
    
    yTemp, urms = plus_invert(yPlus, urms_plus, uTau, kinVisc, Ny)
    yTemp, vrms = plus_invert(yPlus, vrms_plus, uTau, kinVisc, Ny)
    yTemp, wrms = plus_invert(yPlus, wrms_plus, uTau, kinVisc, Ny)
    yTemp, uv = plus_invert(yPlus, uv_plus, uTau, kinVisc, Ny)
    
    RStress = np.zeros((Ny, 4))
    
    for i in range(Ny):
        RStress[i, 0] = urms[i]
        RStress[i, 1] = vrms[i]
        RStress[i, 2] = wrms[i]
        RStress[i, 3] = uv[i]
    
    return RStress    
    
def amplitudes(RStress):
    '''
    Function to calculate the values for the A matrix used in the Lund method
    
    Takes an input of 4 Reynolds stress values, retuns an the full Lund array
    '''
    
    A1 = RStress[0]
    A4 = (abs(RStress[3])/RStress[3])*RStress[3]**2/RStress[0]
    A2 = (RStress[1]**2 - A4**2)**0.5
    A3 = RStress[2]
    
    return np.array([A1, A2, A3, A4])

def blend_amplitudes(inputs):

    # gather global variables
    configuration = inputs.configuration
    Ny = inputs.Ny
    Nz = inputs.Nz
    yGrid = inputs.yGrid
    zGrid = inputs.zGrid

    # Get mean profiles for bottom and sidewalls
    Umean_y, Vmean_y, Tmean_y, Dmean_y = mean_inflow_profile(inputs)
    Umean_z, Vmean, Tmean_z, Dmean_z = mean_inflow_profile(inputs, direction='z')

    # process reynolds stresses
    RStress = process_fluctuation_input(inputs, Umean_y, Tmean_y, Dmean_y)
    RStress_z = process_fluctuation_input(inputs, Umean_z, Tmean_z, Dmean_z, direction='z')

    # initialise amplittudes
    Az = np.zeros((Nz, 6))
    Ay = np.zeros((Ny, 6))

    # boundary layer thickness
    blt = np.interp(0.99, Umean_y, yGrid)

    # get amplitudes from reynolds stress data
    for j in range(Ny):
        for k in range(Nz):
            a = amplitudes(RStress[j])
            Ay[j, :] = [a[0], a[3], a[1], 0, 0, a[2]]
            a = amplitudes(RStress_z[k])
            Az[k, :] = [a[0], 0, a[2], a[3], 0, a[1]]

    if any(conf in configuration for conf in ["full", "parallel"]):
        # full configuration copies both sidewall and bottom wall distribution
        mid_height = int(np.ceil(Ny/2))
        for j in range(mid_height, Ny):
            Ay[j, :] = Ay[Ny - j - 1, :]
            Umean_y[j] = Umean_y[Ny - j - 1]
    # process for copying amplitude distribution on to opposite walls - depends on sidewall configuration
    if any(conf in configuration for conf in ["bottom", "parallel"]):
        # Bottom wall only just copies Ay distribution
        A = np.repeat(Ay[:,np.newaxis,:], Nz, axis=1)
    else:
        if any(conf in configuration for conf in ["full", "channel"]):
            # channel configuration only copies sidewall distribution
            mid_span = int(np.ceil(Nz/2))
            for k in range(mid_span, Nz):
                Az[k, :] = Az[Nz - k - 1, :]
                Umean_z[k] = Umean_z[Nz - k - 1]

        # Blending method
        # ylim, zlim determine the ranges over which the blending method is performed (either halfway or full width)
        if configuration=='corner':
            ylim = int(np.interp(0.50*max(yGrid), yGrid, np.linspace(0, Ny-1, Ny)))
            zlim = int(np.interp(0.50*max(zGrid), zGrid, np.linspace(0, Nz-1, Nz)))
        elif 'channel' in configuration:
            ylim = int(np.interp(0.50*max(yGrid), yGrid, np.linspace(0, Ny-1, Ny)))
            zlim = Nz - 1
        elif configuration=='full':
            ylim = Ny - 1
            zlim = Nz - 1

        # initialise A as a summation of Az and Ay (for the far-field values)
        A = np.repeat(Ay[:,np.newaxis,:], Nz, axis=1) + np.repeat(Az[np.newaxis, :,:], Ny, axis=0)

        # loop within the limits
        for k in range(zlim + 1):
            for j in range(ylim + 1):

                # get the distances to the nearest wall
                ydist = min(yGrid[j], yGrid[-1] - yGrid[j])
                zdist = min(zGrid[k], zGrid[-1] - zGrid[k])

                # corner parameter, h: used to determine weight parameters c1 and c2
                h = ydist - zdist

                # outwith corner:
                if h > blt:
                    c1 = 0.0
                elif h < -blt:
                    c1 = 1.0
                else:
                    # within corner region, compute weight parameters
                    c1 = 0.5 - 0.5*math.sin(0.5*pi*h/blt)

                c2 = 1.0 - c1

                if zdist <= blt:
                    c1 *= (np.absolute(zdist)/blt)**0.25
                if ydist <= blt:
                    c2 *= (np.absolute(ydist)/blt)**0.25

                # Do blending to compute A
                A[j, k, :] = c1*Ay[j, :] + c2*Az[k, :]

    return A

def get_all_vars(constants, values):
    return input_data(constants, values)    
    
        
def inflow_init(constants, values, wall_normal_v=True, fixed_seed=False):
    '''
    Method called within setup file which runs the pre-alrothim steps:
    - stores inputs
    - computes mean profile:
    - processes target profile (amplitude) data
    '''
    inputs = input_data(constants, values)
    inputs.set_wall_normal_v(wall_normal_v)
    Ny = inputs.Ny
    Nz = inputs.Nz
    niter = inputs.Nt
    Umean, Vmean, Wmean, Tmean, Dmean = blend_mean_profile(inputs)
    A = blend_amplitudes(inputs)
        
    seeds = get_seeds(Ny, Nz, niter, fixed_seed=fixed_seed)
    
    return Umean, Vmean, Wmean, Tmean, Dmean, A, seeds
