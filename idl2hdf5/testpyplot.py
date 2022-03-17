import h5py
import numpy as np
import matplotlib.pyplot as plt

filename = "../Data/HDF5/testdata.0010.h5"

with h5py.File(filename, "r") as f:
    # List all groups
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]

    # Get the data
    #data = list(f[a_group_key])
    data=list(f['ro_p'])
    xg=list(f['xgrid'])

nx=len(data)
#ny=len(data[0])

xlist = np.linspace(0.0, 1.0, nx)
#ylist = np.linspace(0.0, 1.0, ny)
#X, Y = np.meshgrid(xlist, ylist)

fig,ax=plt.subplots(1,1)
cp=ax.plot(xg,data)

plt.show()
