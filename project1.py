import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("test.txt", skiprows = 1)

n = data.size

y = np.zeros(n)
x = np.zeros(n)

plt.plot(data)
plt.show()
