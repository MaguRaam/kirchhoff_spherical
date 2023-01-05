import numpy as np 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from math import sqrt

############################################################################################### 

p_inf     = 100.0     # Pressure of liquid far away from the bubble 
rho_inf   = 1.0       # Density of liquid 
gamma_inf = 4.4       # Ratio of specific heats of the liquid 
pi_inf    = 6000.0    # Stiffness constant of the liquid 
c_inf = sqrt(gamma_inf*(p_inf+pi_inf)/rho_inf) # Speed of sound in the liquid 

print("speed of sound in water = ", c_inf)

############################################################################################### 

p0    = 1           # Initial pressure inside the bubble  
gamma = 1.4 		# Ratio of specific heats inside the bubble 
R0    = 1.0         # Initial radius of the bubble 
R0dot = 0.0         # Initial velocity of the bubble interface 

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
    
    term1 = (1.0 - y[1]/c_inf)*y[0]
    term2 = ((3./2.)*(1. - y[1]/(3.0*c_inf)))*y[1]*y[1]
    term3 = ((1.0 + y[1]/c_inf)*(p_B-p_inf) + (y[0]/c_inf)*dp_Bdt)/rho_inf
    
    dydt[1] = (term3-term2)/term1    

    return dydt


t0 = 0.0 
tf = 1.0
y0 = np.array([R0, R0dot])
t = np.linspace(t0,tf,1000)
sol = odeint(kellerMiksis, y0, t)
R = sol[:,0]
Rdot = sol[:,1]

plt.plot(t,R, color='blue', label='Keller-Kolodner')
plt.xlabel(r'$t$')
plt.ylabel(r'$R$')
plt.title('Bubble Collapse')
plt.legend()
plt.savefig('ref.png',dpi=500)

##############################################################################################



