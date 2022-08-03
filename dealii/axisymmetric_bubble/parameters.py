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
R  =  0.25

#Grid: (z, r)
zmin    =  -30.0*R
zmax    =   30.0*R
rmin    =   0
rmax    =   30.0*R
Nz      =   2000
Nr      =   1000
h       =   (zmax - zmin)/Nz
Tf      =   0.15

print("cell size = ", h)
print("no of cells/Radius of the bubble = ", R/h)

#spherical kirchhoff surface
Rk       =   25.0*R                         #radius of semicircle arc
Lk       =   math.pi*Rk                     #length of semicircle arc
Nk       =   5000                           #no of cells on arc

print("no of cells on semicircle arc = ", Nk)

#observer point location:
zo = 0.
ro = 29.0*R


#speed of sound
print("speed of sound in water = ", math.sqrt(g1*(p1 + p_inf1)/rho1) );



