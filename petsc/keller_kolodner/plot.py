import numpy as np
import matplotlib.pyplot as plt

t, R, p = np.loadtxt('euler.dat', unpack=True)
t_euler, p_euler, pder_euler = np.loadtxt('cfd.dat', unpack=True)

plt.plot(t, R, color='blue', label='Keller-Kolodner')
plt.xlabel(r'$t$')
plt.ylabel(r'$R$')
plt.title('Bubble collapse')
plt.legend()
plt.savefig('radius.png', dpi=500)
plt.clf()

plt.plot(t, p, color='blue', label='Keller Kolodner')
plt.plot(t_euler, p_euler, color='red', label='Euler')
plt.xlabel(r'$t$')
plt.ylabel(r'$p - p_{\infty}$')
plt.title('Far-field pressure')
plt.legend()
plt.savefig('pressure.png', dpi=500)
plt.clf()
