import numpy as np
import matplotlib.pyplot as plt


def plot(filename):
    data = np.loadtxt(filename, skiprows = 1)

    x = data[:,0]
    numeric = data[:,1]
    exact = data[:,2]

    n = x.size

    

    plt.xlabel(r"$x$")
    plt.xlabel(r"$x$")
    plt.title("n =" + str(n))
    plt.plot(x, numeric, label ="Numeric")
    plt.plot(x, exact, label = "Exact")
    plt.legend()
    plt.savefig(str(n) + ".png")
    plt.close()

plot("Bn1.txt")
plot("Bn2.txt")
plot("Bn3.txt")

