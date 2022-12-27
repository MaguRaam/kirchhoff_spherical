import numpy as np
import matplotlib.pyplot as plt
import os


def plot(fileName):
    data = np.loadtxt(fileName, delimiter=',')

    x = data[:, 0]
    rho = data[:, 1]
    u = data[:, 2]
    p = data[:, 3]
    phi = data[:, 4]
    p_acoustic = np.subtract(p, 1.0)

		#plot velocity	
	
    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()
    ax1.plot(x, u, 'b-')
    ax2.plot(x, phi, 'r-')

    ax1.set_xlabel('x')
    ax1.set_ylabel('velocity', color='b')
    ax2.set_ylabel(r'$\alpha$', color='r')

    outFile = fileName.replace("csv", "png")
    outFile = outFile.replace("sol", "u")
    plt.savefig(outFile, dpi=400)
    plt.clf()
    
    #plot acoustic pressure:
    
    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()
    ax1.plot(x, p_acoustic, 'b-')
    ax2.plot(x, phi, 'r-')

    ax1.set_xlabel('x')
    ax1.set_ylabel('acoustic pressure $p - p_0$', color='b')
    ax2.set_ylabel(r'$\alpha$', color='r')

    outFile = fileName.replace("csv", "png")
    outFile = outFile.replace("sol", "p")
    plt.savefig(outFile, dpi=400)
    plt.clf()
    


for fileName in os.listdir('.'):
    if fileName.endswith('.csv'):
    	plot(fileName)
