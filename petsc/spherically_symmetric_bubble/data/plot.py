import numpy as np
import matplotlib.pyplot as plt
import os


def plot(fileName):
    data = np.loadtxt(fileName)

    t = data[:, 0]
    p_acoustic = data[:, 1]
    
    plt.plot(t, p_acoustic, 'b-')
    plt.set_xlabel('t')
    plt.set_ylabel('acoustic pressure, $p - p_0$', color='b')

    outFile = fileName.replace("dat", "png")
    plt.savefig(outFile, dpi=400)
    plt.clf()


for fileName in os.listdir('.'):
    print(fileName)
    plot(fileName)
