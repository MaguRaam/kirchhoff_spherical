#simulation parameters!
from numpy import double
import math

#----------------------Scaled Values----------------------------------------#

#properties of multiphase system (air bubble in water)

#water medium  (phi = 1):
g1                = 4.4          #Specific heat ratio of water
p_inf1            = 6.0e8        #Stiffness constant of water
rho1              = 1000.0       #initial density of water
p1                = 101325.0     #initial pressure of water


#air bubble (phi = 0):
g2                = 1.4          #Specific heat ratio of air
p_inf2            = 0.0          #Stiffness constant of air
rho2              = 1.2246       #initial density of air
p2                = 101325.0     #initial pressure of air

#------------------------------------------------------#

#bubble radius:
R  =  6.0e-5

#Grid: (z, r)
zmin    =  -6.6e-4
zmax    =   3.4e-4
rmin    =   0
rmax    =   5.e-4
Nz      =   392
Nr      =   196
h       =   (zmax - zmin)/Nz

print("cell size = ", h)
print("no of cells/Radius of the bubble = ", R/h)

#spherical kirchhoff surface
Rk       =   5.0*R                           #radius of semicircle arc
Lk       =   math.pi*Rk                      #length of semicircle arc
Nk       =   90                              #no of quadrature points

#speed of sound
print("speed of sound in water = ", math.sqrt(g1*(p1 + p_inf1)/rho1) );


