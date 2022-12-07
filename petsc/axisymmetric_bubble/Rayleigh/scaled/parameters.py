#simulation parameters!
from numpy import double
import math

#----------------------Scaled Values----------------------------------------#

#properties of multiphase system (air bubble in water)

#water medium  (phi = 1):
g1                = 4.4          #Specific heat ratio of water
p_inf1            = 6000.0       #Stiffness constant of water
rho1              = 1.0       	 #initial density of water  (scaled by 1e3) so 1000kgm-3
p1                = 10.0         #initial pressure of water (scaled by 1e5) so 1MPa


#air bubble (phi = 0):
g2                = 1.4          #Specific heat ratio of air
p_inf2            = 0.0          #Stiffness constant of air
rho2              = 0.001        #initial density of air (scaled by 1e3) so 1Kgm-3
p2                = 1.0          #initial pressure of air (scaled by 1e5) so .1MPa

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
Nk       =   500

print("no of cells on semicircle arc = ", Nk)

#observer point location:
io = 250
jo = 225

print("zo = ", zmin + (io + 0.5)*h)
print("ro = ", rmin + (jo + 0.5)*h)
print("Kirchhoff surface Rk = ", Rk)


#speed of sound
print("speed of sound in water = ", math.sqrt(g1*(p1 + p_inf1)/rho1) );


