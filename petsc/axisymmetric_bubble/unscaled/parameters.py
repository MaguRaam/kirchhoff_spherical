#simulation parameters!
from numpy import double
import math

#----------------------Scaled Values----------------------------------------#

#properties of multiphase system (air bubble in water)

#water medium  (phi = 1):
g1                = 4.4          #Specific heat ratio of water
p_inf1            = 6000.0       #Stiffness constant of water
rho1              = 1000.0       #initial density of water
p1                = 1.0e6        #initial pressure of water


#air bubble (phi = 0):
g2                = 1.4          #Specific heat ratio of air
p_inf2            = 0.0          #Stiffness constant of air
rho2              = 1.0          #initial density of air
p2                = 1.0e5        #initial pressure of air

#------------------------------------------------------#

#bubble radius:
R  =  0.038

#Grid: (z, r)
zmin    =  -10.0*R
zmax    =   10.0*R
rmin    =   0
rmax    =   10.0*R
Nz      =   500
Nr      =   250
h       =   (zmax - zmin)/Nz

print("cell size = ", h)
print("no of cells/Radius of the bubble = ", R/h)

#spherical kirchhoff surface
Rk       =   6.0*R                           #radius of semicircle arc
Lk       =   math.pi*Rk                     #length of semicircle arc


#observer point location:
io = 250
jo = 225

print("zo = ", zmin + (io + 0.5)*h)
print("ro = ", rmin + (jo + 0.5)*h)
print("Kirchhoff surface Rk = ", Rk)


#speed of sound
print("speed of sound in water = ", math.sqrt(g1*(p1 + p_inf1)/rho1) );


