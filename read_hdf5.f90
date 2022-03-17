PROGRAM read_hdf5
  USE HDF5 ! This module contains all necessary modules
  IMPLICIT NONE

  character(40) :: file_path
  integer(HID_T) :: dset_id, file_id, filespace_id, memspace_id, plist_id
  integer(HSIZE_T) :: dimsFile(1), dimsMem(1)
  integer(HSSIZE_T) :: offsetFile(1)
  integer :: hdferr, request, margin, mpi_size, mpi_pos
  double precision :: data_out(1002)

  request = 1000
  mpi_size = 8
  margin = 1
  mpi_pos = 3
  dimsMem(1) = request/mpi_size + 2*margin
  dimsFile(1) = request + 2*margin
  offsetFile(1) = request*mpi_pos/mpi_size

  file_path = "Data/t0004.h5"
  CALL h5open_f(hdferr)
  CALL h5fopen_f(trim(file_path), H5F_ACC_RDONLY_F, file_id, hdferr)
  !  Data transfer property (default is H5FD_MPIO_INDEPENDENT_F)
  !CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdferr)
  !  Get dataset from file
  CALL h5dopen_f(file_id, "xgrid", dset_id, hdferr)
  !  Get file dataspace
  CALL h5dget_space_f(dset_id, filespace_id, hdferr)
  CALL h5sselect_hyperslab_f(filespace_id, H5S_SELECT_SET_F, offsetFile, dimsMem, hdferr)
  !  Set memory dataspace
  CALL h5screate_simple_f(1, dimsMem, memspace_id, hdferr)
  !CALL h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, offset, dims_1D, hdferr)
  ! read file
  !CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdferr)
  CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_out, dimsFile, hdferr, &
                 mem_space_id=memspace_id, file_space_id=filespace_id)
  ! close connections
  !CALL h5pclose_f(plist_id, hdferr)
  CALL h5sclose_f(filespace_id, hdferr)
  CALL h5sclose_f(memspace_id, hdferr)
  CALL h5dclose_f(dset_id, hdferr)
  CALL h5fclose_f(file_id, hdferr)
  CALL h5close_f(hdferr)

  write(*,*) file_path
  write(*,*) data_out(1), data_out(dimsMem(1))
END PROGRAM read_hdf5
