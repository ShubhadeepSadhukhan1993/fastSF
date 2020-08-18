import numpy as np
import h5py

def hdf5_reader(filename,dataset):
    file_V1_read = h5py.File(filename)
    dataset_V1_read = file_V1_read["/"+dataset]
    V1=dataset_V1_read[:,:,:]
    return V1


U = hdf5_reader("U.V2r.h5", "U.V2r")

f1 = h5py.File("Y.F.h5", "w")
dset1 = f1.create_dataset("Y.F", data = U)
f1.close()
