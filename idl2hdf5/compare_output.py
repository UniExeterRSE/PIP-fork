import h5py
import numpy as np
import pdb

# Looping over time steps
for i in np.arange(11):
    fn_ref = f"../Data/run4/t{format(i, '04d')}.h5"
    #fn_ref = f"testdata_3D_OZ.{format(i, '04d')}.h5"
    # define hdf5 filename of file to check
    fn_new = f"../Data/run3/t{format(i, '04d')}.h5"

    print(f"time step #{i}")
    # For each reference file (contains multiple variables)
    with h5py.File(fn_ref, "r") as ref:
        with h5py.File(fn_new, "r") as new:
            for param in ref:
                ref_data = np.array(ref[param])
                new_data = np.array(new[param])

                # Compare element by element
                comp = np.isclose(ref_data, new_data)
                if not comp.all():
                    print(f"\t{param} parameter array has {comp.size - comp.sum()} mismatches")
                    pdb.set_trace()
