#Plot exact solution and Kirchoff solution at observer location 

import numpy as np
import matplotlib.pyplot as plt


p = np.loadtxt("pressure.dat")


fig,ax = plt.subplots()
plt.figure(figsize=(14,14))
plt.rcParams.update({'font.size': 30})
plt.xlabel('observer time, $t_{0}$')
plt.ylabel('acoustic pressure, $p - p_{0}$')
leg1,= plt.plot(p[:,0], p[:,1],"b",marker = 'o', markersize=10, linewidth=3)
leg2,= plt.plot(p[:,0],p[:,2],"r-", linewidth=3)
plt.legend((leg1, leg2), ('Kirchhoff sol.', 'Analytical sol.'), loc='upper right', fontsize=22, numpoints=1)
plt.savefig("Pressure")