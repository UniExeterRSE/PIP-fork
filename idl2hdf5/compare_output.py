import h5py
import numpy as np
import pdb

# Looping over time steps
for i in np.arange(11):
    fn_dac = "testdata_con.{0:04d}.h5".format(i)
    # define hdf5 filename of file to check
    fn_h5 = f"../Data/t{format(i, '04d')}.h5"

    # For each reference file (contains multiple variables)
    with h5py.File(fn_dac, "r") as ref:
        with h5py.File(fn_h5, "r") as new:
            for param in ref:
                ref_data = np.array(ref[param])
                new_data = np.array(new[param]).flatten()

                # Compare element by element
                comp = np.isclose(ref_data, new_data)
                if not comp.all():
                    print(f"{param} parameter array in time step #{i} has mismatches")
                    #pdb.set_trace()
