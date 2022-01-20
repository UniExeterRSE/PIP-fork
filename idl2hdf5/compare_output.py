import h5py
import numpy as np
import pdb

# Looping over time steps
for i in np.arange(11):
    fn_dac = "testdata_con.{0:04d}.h5".format(i)

    # For each reference file (contains multiple variables)
    with h5py.File(fn_dac, "r") as ref:
        for param in ref:
            ref_data = np.array(ref[param])

            # define hdf5 filename of file to check
            if param == 'xgrid':
                fn_h5 = "../Data/x.dac.0000.h5"
            else:
                fn_h5 = f"../Data/{format(i, '04d')}{param}.dac.0000.h5"

            # open relevant parameter simulaion hdf5 file
            with h5py.File(fn_h5, "r") as new:
                # there is only one key in each current sim file
                only_key = list(new.keys())[0]

                new_data = np.array(new[only_key]).flatten()
                # Compare element by element
                comp = np.isclose(ref_data, new_data)
                if not comp.all():
                    print(f"{param} parameter array in time step #{i} has mismatches")
                    #pdb.set_trace()
