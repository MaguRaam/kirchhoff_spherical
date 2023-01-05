import numpy as np 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from math import sqrt

############################################################################################### 

p_inf     = 100.0     # Pressure of liquid far away from the bubble 
rho_inf   = 1.0     # Density of liquid 
gamma_inf = 4.4     # Ratio of specific heats of the liquid 
pi_inf    = 6000.0   # Stiffness constant of the liquid 
c_inf = sqrt(gamma_inf*(p_inf+pi_inf)/rho_inf) # Speed of sound in the liquid 
mu = 0.0            # Viscosity of the liquid
S = 0.0             # interface surface tension

print("speed of sound in water = ", c_inf)

############################################################################################### 

p0    = 1.0         # Initial pressure inside the bubble  
gamma = 1.4 		# Ratio of specific heats inside the bubble 
R0    = 1.0         # Initial radius of the bubble 
R0dot = 0.0         # Initial velocity of the bubble interface 

Tc = 0.915*R0*sqrt(rho_inf/(p_inf-p0))  # Rayleigh collapse time 

############################################################################################### 

def kellerMiksis(y, t):
    # t    -> time  
    # y[0] -> radius of bubble (R)
    # y[1] -> velocity of bubble interface   (dR/dt)  
    
    dydt = np.zeros(2)
    dydt[0] = y[1]
    
    # Bubble pressure 

    p_B = (p0*(R0/y[0]))**(3.0*gamma) # Pressure inside the bubble
    dp_Bdt = -3.0*gamma*p0*(R0**(3.0*gamma))*(pow(y[0],(-3.0*gamma-1.0)))*y[1]
    p = p_B - 2.0*S/y[0] - 4.0*mu*y[1]/y[0]
    dpdt = dp_Bdt + (2.0*S*y[1] + 4.0*mu*y[1]*y[1])/(y[0]*y[0]);

    term1 = (1.0 - y[1]/c_inf)*y[0] + 4.0*mu/(rho_inf*c_inf);
    term2 = ((3./2.)*(1. - y[1]/(3.0*c_inf)))*y[1]*y[1];
    term3 = ((1.0 + y[1]/c_inf)*(p-p_inf) + (y[0]/c_inf)*dpdt)/rho_inf;
    
    dydt[1] = (term3-term2)/term1    

    return dydt


def rayleighPlesset(y, t):
    # Special case of Keller-Miksis (c_inf -> infinity)
    # t    -> time  
    # y[0] -> radius of bubble (R)
    # y[1] -> velocity of bubble interface   (dR/dt)  
    
    dydt = np.zeros(2)
    dydt[0] = y[1]
    
    # Bubble pressure 

    pB    = p0*(R0/y[0])**(3.0*gamma)

    term1 = (1.0)*y[0]    
    term2 = (3./2.)*y[1]*y[1]
    term3 = (pB-p_inf)/rho_inf
    
    dydt[1] = (term3-term2)/term1    

    return dydt
############################################################################################### 

t0 = 0.0 
tf = 10*Tc
y0 = np.array([R0, R0dot])
t = np.linspace(t0,tf,1000)
sol = odeint(kellerMiksis, y0, t)
R = sol[:,0]
Rdot = sol[:,1]

print("final time tf = ", tf)

tN,RN = np.loadtxt('bubble_radius.csv', delimiter=',', unpack=True)

plt.plot(t/Tc,R/R0, color='black', label='Keller-Miksis')
plt.plot(tN/Tc,RN/R0, marker='o', color='blue', linestyle='None', markersize=3, markerfacecolor='None', label='Multiphase solver')
#plt.grid()
plt.xlabel(r'$t/T_c$')
plt.ylabel(r'$R/R_0$')
plt.title('Bubble Collapse')
plt.legend()
plt.savefig('fig.png',dpi=500)

##############################################################################################



