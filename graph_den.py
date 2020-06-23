import sys
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib

#-------------------------------------------data plot-------------------------------------------------
data = np.loadtxt('density_0000', skiprows=0)
x = data[:,0]
y = data[:,1]
z = data[:,2]
n = (len(data))
#------------------------------------------plot setting----------------------------------------------- 
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plot_max = sys.argv[1]

count = 0
while (count < plot_max):
#  plot
   
   fname = ('density_%04d')%(count)
   data = np.loadtxt(fname, skiprows=0)
   x = data[:,0]
   y = data[:,1]
   z = data[:,2]
   ax.scatter(x, y, z, c='b', marker='o')
   
   figname = ('den_fig_%04d.png')%(count)
   plt.savefig(fname, dpi=840, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None, metadata=None)

   plt.cla()
   count += 100

