import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

filename = "data.txt"
graph = "graph"

data = np.loadtxt(filename, skiprows = 1)
f = open("data.txt", "r")
line = f.readlines(1)

N = int(line[0].strip().split()[0])
M = int(line[0].strip().split()[1])

x = data[:,0]
y = data[:,1]
z = data[:,2]
x = x.reshape((N,M))
y = y.reshape((N,M))
z = z.reshape((N,M))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, z, cmap='inferno')

ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u')
plt.show()