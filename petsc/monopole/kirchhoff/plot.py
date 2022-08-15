#Plot exact solution and Kirchoff solution at observer location 
import numpy as np
import matplotlib.pyplot as plt


#Plot exact solution and Kirchoff solution at observer location 
p_exact = np.loadtxt("p_exact.dat")
p_kirchhoff = np.loadtxt("p_kirchhoff.dat")


fig,ax = plt.subplots()
plt.figure(figsize=(8,8))
plt.rcParams.update({'font.size': 11})
params = {'legend.fontsize': 60,
          'legend.handlelength': 3}
plt.rcParams.update(params)
#plt.title('perturbed pressure as a function of time at observer point')
plt.xlabel('$observer time, t_{0}$')
plt.ylabel('acoustic pressure, $p - p_{0}$')
leg1,= plt.plot(p_kirchhoff[:,0], p_kirchhoff[:,1],'b-',markersize=1)
leg2,= plt.plot(p_exact[:,0],p_exact[:,1],'r--',markersize=1)
plt.legend((leg1, leg2), ('Kirchhoff', 'exact'), loc='upper right', fontsize=10, numpoints=1)
plt.savefig("Pressure")
plt.grid()
