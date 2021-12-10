; IDL routine to save HDF5 files 

;Input directory
fname='../Data'

;Define save directory MUST ALREADY EXIST
savdir=fname+'/HDF5/'

;Variables to save
outvar=['xgrid','ro_p','vx_p','vy_p','vz_p','bx','by','bz','pr_p']

;Read in the data
for tread=0,30 do begin
;tread=30
rdmpi,ds,datapath=fname,time_step=tread

;Define simulation grid
x=ds.x
ro=ds.ro_p
vx=ds.vx_p
vy=ds.vy_p
vz=ds.vz_p
bx=ds.bx
by=ds.by
bz=ds.bz
pr_p=ds.pr_p

    ;; get data type and space, needed to create the dataset

;    file = savdir+outvar(i)+'.h5'
    file = savdir+'testdata.'+string(tread,FORMAT='(I4.4)')+'.h5'
    fid = H5F_CREATE(file)

    ;datatype_id = H5T_IDL_CREATE(ro)
    ;dataspace_id = H5S_CREATE_SIMPLE(size(ro_p,/DIMENSIONS))

for i=0,n_elements(outvar)-1 do begin
;for i=0,1 do begin

    if outvar(i) eq 'xgrid' then data = x
    if outvar(i) eq 'ro_p' then data = ro
    if outvar(i) eq 'vx_p' then data = vx
    if outvar(i) eq 'vy_p' then data = vy
    if outvar(i) eq 'vz_p' then data = vz
    if outvar(i) eq 'bx' then data = bx
    if outvar(i) eq 'by' then data = by
    if outvar(i) eq 'bz' then data = bz
    if outvar(i) eq 'pr_p' then data = pr_p

    datatype_id = H5T_IDL_CREATE(data)
    dataspace_id = H5S_CREATE_SIMPLE(size(data,/DIMENSIONS))
    
    dataset_id = H5D_CREATE(fid,$

    outvar(i),datatype_id,dataspace_id)

    ;; write data to dataset

    H5D_WRITE,dataset_id,data

    ;; close all open identifiers

    H5D_CLOSE,dataset_id

    H5S_CLOSE,dataspace_id

    H5T_CLOSE,datatype_id

endfor

H5F_CLOSE,fid

endfor

END
