import numpy as np
import matplotlib.pyplot as plt

def error():
    data = []
    x = []

    for i in range(1, 8):
        d = open("Cn" + str(i) + ".txt")
        line = d.readline()
        data.append(float(line))
        h = -i
        x.append(h)
    return data, x

def plot(filename):
    data = np.loadtxt(filename, skiprows = 2)

    x = data[:,0]
    numeric = data[:,1]
    exact = data[:,2]

    n = x.size

    

    plt.xlabel(r"$x$")
    #plt.ylabel(r"$x$")
    plt.title("n =" + str(n))
    plt.plot(x, numeric, "r--",label ="Numeric", )
    plt.plot(x, exact, label = "Exact")
    plt.legend()
    plt.savefig(str(n) + ".png")
    plt.close()

data, x = error()

print(data, x)

plt.xlabel(r"$log_{10}(h)$")
#plt.ylabel(r"$x$")
plt.title("Max relative error ")
plt.plot(x, data, "ro")
plt.legend()
#plt.show()
plt.savefig("ERROR.png")

#plot("Cn1.txt")
#plot("Cn2.txt")
#plot("Cn3.txt")

