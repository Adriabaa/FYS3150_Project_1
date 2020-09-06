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
    plt.show()

plot("big.txt")
plot("bigger.txt")
plot("biggest.txt")

