!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine WriteHDF5(fname,save_pressure)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use hdf5

include 'header'

character*55 fname
logical final,save_pressure

real*8 tmp(nx,ny_s+1,nz_s+1)
real*8 tmp_th(nx_th,ny_s_th+1,nz_s_th+1)
!
!     HDF5 ------------------------------------------------------
!
!     Dataset names
character(len=10) :: dname

!     Identifiers
integer(hid_t) :: file_id, dset_id
integer(hid_t) :: filspace_id, memspace_id
integer(hid_t) :: filspace_id_th, memspace_id_th

!     Identifiers
integer(hid_t) :: selspace_id
integer(hid_t) :: plist_id_w,plist_id_d,plist_id_d_th

!     Dimensions in the memory and in the file
integer(hsize_t), dimension(3) :: dimsm,dimsf
integer(hsize_t), dimension(3) :: dimsm_th,dimsf_th

integer(hsize_t), dimension(3) :: chunk_dims, count, offset_f
integer(hsize_t), dimension(3) :: stride, block, offset_m
integer(hsize_t), dimension(3) :: chunk_dims_th, offset_th
integer(hsize_t), dimension(3) :: block_th

integer :: rHDF5 = 3, arank, i, j, k

integer(hsize_t),dimension(2) :: adims
integer(size_t)               :: tdim
integer(hid_t)                :: aid,tspace,ttype
integer, dimension(1)         :: tint(6)
character*80                  :: namnbuf
character*20                  :: sttimec

integer error, ith

double precision En(4)

dimsm(1:3) = (/nx,ny_s+1,nz_s+1/)
dimsf(1:3) = (/nx,ny,nz/)
block = dimsm
dimsm_th(1:3) = (/nx_th,ny_s_th+1,nz_s_th+1/)
dimsf_th(1:3) = (/nx_th,ny_th,nz_th/)
block_th = dimsm_th

if (rank == 0) write(*,*) 'Writing flow to ',fname
!write(*,*) 'ranky,rankz: ',ranky,rankz

chunk_dims = (/nx,1,nz_s+1/)
chunk_dims_th = (/nx_th,1,nz_s_th+1/)

! Stride and count for number of rows and columns in each dimension
stride = 1
count  = 1

! Offset determined by the rank of a processor
offset_f = (/ 0, ranky*(ny_s+1), rankz*(nz_s+1) /)
offset_th = (/ 0, ranky*(ny_s_th+1), rankz*(nz_s_th+1) /)
offset_m(1:3)=0

! Initialize interface
call h5open_f(error)

! Setup file access property list with parallel I/O access
call h5pcreate_f(h5p_file_access_f, plist_id_d, error)
call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_world, &
                            mpi_info_null, error)

! Create the file collectively
call h5fcreate_f(trim(fname), h5f_acc_trunc_f, &
                    file_id, error, access_prp = plist_id_d)
call h5pclose_f(plist_id_d, error)

adims=(/3,2/)
arank=2
call h5screate_simple_f(arank,adims,tspace, error)

! -----------------------------
! Resolution
call h5acreate_f(file_id,'Resolution',h5t_std_i32le,tspace, &
                        aid, error)
tint(1)=nx
tint(2)=ny
tint(3)=nz
tint(4)=nx_th
tint(5)=ny_th
tint(6)=nz_th
call h5awrite_f(aid,h5t_native_integer,tint,adims,error)
call h5aclose_f(aid, error)
! -----------------------------
! Close dataspace
call h5sclose_f(tspace, error)

! -----------------------------
! Date

adims=(/1,20/)
arank=1
tdim=20
call h5tcopy_f(h5t_fortran_s1, ttype, error)
call h5tset_size_f(ttype,tdim,error)
call h5screate_simple_f(arank,adims,tspace, error)

call h5acreate_f(file_id,'Date',ttype,tspace,aid, error)
call time_string(sttimec)
call h5awrite_f(aid,ttype,trim(sttimec),adims, error)
call h5aclose_f(aid, error)
call h5sclose_f(tspace,error)
call h5tclose_f(ttype,error)
! -----------------------------

! -----------------------------
! Extra info
adims=(/1,80/)
tdim=80
call h5tcopy_f(h5t_fortran_s1,ttype,error)
call h5tset_size_f(ttype,tdim,error)
call h5screate_simple_f(arank,adims,tspace, error)

call h5acreate_f(file_id,'Info',ttype,tspace,aid, error)
namnbuf=' (Put here what you want) '
call h5awrite_f(aid,ttype,trim(namnbuf),adims, error)
call h5aclose_f(aid, error)
call h5sclose_f(tspace,error)
call h5tclose_f(ttype,error)
! -----------------------------

call h5screate_f(h5s_scalar_f,tspace,error)

! -----------------------------
! Time & Timestep
call h5acreate_f(file_id,'Time',h5t_ieee_f64le,tspace, &
                        aid, error)
call h5awrite_f(aid,h5t_native_double,time,adims,error)
call h5aclose_f(aid, error)

call h5acreate_f(file_id,'Timestep',h5t_std_i32le,tspace, &
                        aid, error)
call h5awrite_f(aid,h5t_native_integer,time_step,adims,error)
call h5aclose_f(aid, error)
! ----------------------------

call h5sclose_f(tspace,error)

! Convert to physical space
call fft_xzy_mpi_to_physical(cu1,u1)
call fft_xzy_mpi_to_physical(cu2,u2)
call fft_xzy_mpi_to_physical(cu3,u3)
do ith=1,n_th
    csth1(:,:,:)=cth(:,:,:,ith)
    call fft_xzy_mpi_to_physical_th(csth1,sth1)
    th(:,:,:,ith)=sth1(:,:,:)
end do

! Create property list for the chunked dataset creation
call h5pcreate_f(h5p_dataset_create_f, plist_id_d, error)
call h5pset_chunk_f(plist_id_d, rHDF5, chunk_dims, error)
call h5pcreate_f(h5p_dataset_create_f, plist_id_d_th, error)
call h5pset_chunk_f(plist_id_d_th, rHDF5, chunk_dims_th, error)

! Create the dataspaces for velocities and scalars
call h5screate_simple_f(rHDF5, dimsf, filspace_id, error)
call h5screate_simple_f(rHDF5, dimsm, memspace_id, error)
call h5screate_simple_f(rHDF5, dimsf_th, filspace_id_th, error)
call h5screate_simple_f(rHDF5, dimsm_th, memspace_id_th, error)

! Create property list for collective dataset write
call h5pcreate_f(h5p_dataset_xfer_f, plist_id_w, error)
call h5pset_dxpl_mpio_f(plist_id_w, h5fd_mpio_collective_f, error)

do ith=1,3+n_th

    select case(ith)
    case (1)
        do j=0,ny_s
            do k=0,nz_s
                do i=0,nxm
                    tmp(i+1,j+1,k+1)=u1(i,k,j)
                end do
            end do
        end do
        dname="u1"
    case (2)
        do j=0,ny_s
            do k=0,nz_s
                do i=0,nxm
                    tmp(i+1,j+1,k+1)=u2(i,k,j)
                end do
            end do
        end do
        dname="u2"
    case (3)
        do j=0,ny_s
            do k=0,nz_s
                do i=0,nxm
                    tmp(i+1,j+1,k+1)=u3(i,k,j)
                end do
            end do
        end do
        dname="u3"
    case (4:)
        do j=0,ny_s_th
            do k=0,nz_s_th
                do i=0,nxm_th
                    tmp_th(i+1,j+1,k+1)=th(i,k,j,ith-3)
                end do
            end do
        end do
        dname="th"//char(ith+45)
    end select
    if (ith <= 3) then
! If saving velocity, use velocity grid
        call h5dcreate_f(file_id, trim(dname), h5t_ieee_f64le, &
                filspace_id, dset_id, error, dcpl_id = plist_id_d)

! Select hyperslab in the file.
!           call h5dget_space_f(dsetur_id, selspace_id, error)
            call h5sselect_hyperslab_f (filspace_id, h5s_select_set_f, &
                    offset_f, count, error, stride, block)

            call h5sselect_hyperslab_f (memspace_id, h5s_select_set_f, &
                    offset_m, count, error, stride, block)

! Write the dataset collectively
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                tmp, &
                dimsm, error, file_space_id = filspace_id, &
                mem_space_id = memspace_id, xfer_prp = plist_id_w)

! Close dateset
        call h5dclose_f(dset_id, error)
    else
! If saving scalar, use TH grid
        call h5dcreate_f(file_id, trim(dname), H5T_IEEE_F64LE, &
                filspace_id_th, dset_id, error, dcpl_id = plist_id_d_th)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
            call h5sselect_hyperslab_f (filspace_id_th, H5S_SELECT_SET_F, &
                    offset_th, count, error, stride, block_th)

            call h5sselect_hyperslab_f (memspace_id_th, H5S_SELECT_SET_F, &
                    offset_m, count, error, stride, block_th)

! Write the dataset collectively
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                    tmp_th, &
                    dimsm_th, error, file_space_id = filspace_id_th, &
                    mem_space_id = memspace_id_th, xfer_prp = plist_id_w)

! Close dateset
        call h5dclose_f(dset_id, error)
    end if
end do

!     In the case of saving for the pressure as well
if (save_pressure) then
    call fft_xzy_mpi_to_physical(cp,p)

    do j=0,ny_s
        do k=0,nz_s
            do i=0,nxm
                tmp(i+1,j+1,k+1)=p(i,k,j)
            end do
        end do
    end do
    dname="p"

    call h5dcreate_f(file_id, trim(dname), h5t_ieee_f64le, &
                filspace_id, dset_id, error, dcpl_id = plist_id_d)

!     Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
        call h5sselect_hyperslab_f (filspace_id, h5s_select_set_f, &
                    offset_f, count, error, stride, block)

        call h5sselect_hyperslab_f (memspace_id, h5s_select_set_f, &
                    offset_m, count, error, stride, block)

!     Write the dataset collectively
        call h5dwrite_f(dset_id, h5t_native_double, &
                    tmp, &
                    dimsm, error, file_space_id = filspace_id, &
                    mem_space_id = memspace_id, xfer_prp = plist_id_w)

!     Close dateset
    call h5dclose_f(dset_id, error)
    call fft_xzy_mpi_to_fourier(p,cp)
end if

!     Close the dataspace for the memory and for the file
call h5sclose_f(filspace_id, error)
call h5sclose_f(memspace_id, error)
call h5sclose_f(filspace_id_th, error)
call h5sclose_f(memspace_id_th, error)

!     Close the properties for the dataspace creation and the writing
call h5pclose_f(plist_id_d, error)
call h5pclose_f(plist_id_d_th, error)
call h5pclose_f(plist_id_w, error)

      ! Close groups
call h5fclose_f(file_id, error)
call h5close_f(error)

call fft_xzy_mpi_to_fourier(u1,cu1)
call fft_xzy_mpi_to_fourier(u2,cu2)
call fft_xzy_mpi_to_fourier(u3,cu3)
do ith=1,n_th
    sth1(:,:,:)=th(:,:,:,ith)
    call fft_xzy_mpi_to_fourier_th(sth1,csth1)
    cth(:,:,:,ith)=csth1(:,:,:)
end do

! call mpi_finalize(ierror)
! stop

end subroutine WriteHDF5

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine time_string(cdt)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

!     Construct string in the format '19-DEC-2005 22:47:06'

implicit none

integer i

integer val(8)
character*20 cdt
character*3 monc
character*2 day, hour, minute, sec
character*4 year

call date_and_time(values=val)

if (val(2) == 1) then
    monc  = 'JAN'
else if (val(2) == 2) then
    monc  = 'FEB'
else if (val(2) == 3) then
    monc  = 'MAR'
else if (val(2) == 4) then
    monc  = 'APR'
else if (val(2) == 5) then
    monc  = 'MAY'
else if (val(2) == 6) then
    monc  = 'JUN'
else if (val(2) == 7) then
    monc  = 'JUL'
else if (val(2) == 8) then
    monc  = 'AUG'
else if (val(2) == 9) then
    monc  = 'SEP'
else if (val(2) == 10) then
    monc  = 'OCT'
else if (val(2) == 11) then
    monc  = 'NOV'
else if (val(2) == 12) then
    monc  = 'DEC'
else
    monc  = 'XXX'
end if

write(day,'(i0.2)') val(3)
write(year,'(i4)') val(1)
write(hour,'(i0.2)') val(5)
write(minute,'(i0.2)') val(6)
write(sec,'(i0.2)') val(7)

cdt = day // '-' // monc // '-' // year // ' ' // &
        hour // ':' // minute // ':' // sec

end subroutine time_string

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine ReadHDF5(fname)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use hdf5

include 'header'

character*55 fname
logical final, read_pressure

real*8 tmp(nx,ny_s+1,nz_s+1)
real*8 tmp_th(nx_th,ny_s_th+1,nz_s_th+1)

!     Dataset names
character(len=10) :: dname

!     Identifiers
integer(hid_t) :: file_id, dset_id
integer(hid_t) :: filspace_id, memspace_id
integer(hid_t) :: filspace_id_th, memspace_id_th

!     Identifiers
integer(hid_t) :: selspace_id
integer(hid_t) :: plist_id_w,plist_id_d

!     Dimensions in the memory and in the file
integer(hsize_t), dimension(3) :: dimsm
integer(hsize_t), dimension(3) :: dimsm_th

integer(hsize_t), dimension(3) :: count, offset_f
integer(hsize_t), dimension(3) :: stride, block, offset_m
integer(hsize_t), dimension(3) :: offset_th, block_th

integer :: rHDF5 = 3

integer(hsize_t),dimension(2) :: adims
integer(hid_t)                :: aid
integer, dimension(2)         :: tint(3,2)
character*80                  :: namnbuf
character*20                  :: sttimec

integer error, ith, i, j, k

double precision En(4)

dimsm(1:3) = (/nx,ny_s+1,nz_s+1/)
block = dimsm
dimsm_th(1:3) = (/nx_th,ny_s_th+1,nz_s_th+1/)
block_th = dimsm_th

!     Stride and count for number of rows and columns in each dimension
stride = 1
count  = 1

!     Offset determined by the rank of a processor
offset_f = (/ 0, ranky*(ny_s+1), rankz*(nz_s+1) /)
offset_th = (/ 0, ranky*(ny_s_th+1), rankz*(nz_s_th+1) /)
offset_m(1:3)=0

!     Initialize interface
call h5open_f(error)

!     Setup file access property list with parallel I/O access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id_d, error)
call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_world, &
                        mpi_info_null, error)

!     Open the file
call h5fopen_f(trim(FNAME), H5F_ACC_RDONLY_F, &
                file_id, error, access_prp = plist_id_d)
call h5pclose_f(plist_id_d, error)

!     -----------------------------
!     Resolution
!     -----------------------------
adims=(/3,2/)
call h5aopen_by_name_f(file_id,'.','resolution',aid,error)
call h5aread_f(aid,h5t_native_integer,tint,adims,error)
call h5aclose_f(aid, error)
! Check that the resolution is of the same kind
if ((tint(1,1) /=  nx)  .or.  (tint(2,1) /= ny)  .or. &
        (tint(3,1) /= nz) .or. (tint(1,2) /= nx_th) .or. &
        (tint(2,2) /= ny_th) .or. (tint(3,2) /= nz_th)) then
    if (rank == 0) then
        write(*,*) ' Error. File and program have &
                    &different resolutions. '
        write(*,*) ' Program: ', nx,ny,nz,nx_th,ny_th,nz_th
        write(*,*) ' File   : ', tint(1:3,1),tint(1:3,2)
    end if
    call mpi_finalize(ierror)
    stop
end if

!     -----------------------------
!     Time stamps
!     -----------------------------
adims=(/1,80/)
call h5aopen_by_name_f(file_id,'.','Time',aid,error)
call h5aread_f(aid,h5t_native_double,time,adims,error)
call h5aclose_f(aid, error)

call h5aopen_by_name_f(file_id,'.','Timestep',aid,error)
call h5aread_f(aid,h5t_native_integer,time_step,adims,error)
call h5aclose_f(aid, error)
!     -----------------------------

!     Create property list for collective dataset read
call h5pcreate_f(H5P_DATASET_XFER_F, plist_id_w, error)
call h5pset_dxpl_mpio_f(plist_id_w, H5FD_MPIO_COLLECTIVE_F, error)
!     Dataspaces in memory
call h5screate_simple_f(rHDF5, dimsm, memspace_id, error)
call h5screate_simple_f(rHDF5, dimsm_th, memspace_id_th, error)

do ith=1,3+n_th
!     Here it starts the loop--->
    select case(ith)
    case (1)
        dname="u1"
    case (2)
        dname="u2"
    case (3)
        dname="u3"
    case (4:)
        dname="th"//char(ith+45)
    end select

    if (ith <= 3) then
        call h5dopen_f(file_id,trim(dname),dset_id,error)
        call h5dget_space_f(dset_id,filspace_id,error)
! Select hyperslab in the file.
        call h5sselect_hyperslab_f (filspace_id, h5s_select_set_f, &
                        offset_f, count, error, stride, block)
        call h5sselect_hyperslab_f (memspace_id, h5s_select_set_f, &
                        offset_m, count, error, stride, block)
! Read the dataset collectively
        call h5dread_f(dset_id, h5t_native_double, tmp, &
                    dimsm, error, file_space_id = filspace_id, &
                    mem_space_id = memspace_id, xfer_prp = plist_id_w)
        select case(ith)
        case (1)
            do j=0,ny_s
                do k=0,nz_s
                    do i=0,nxm
                        u1(i,k,j)=tmp(i+1,j+1,k+1)
                    end do
                end do
            end do
        case (2)
            do j=0,ny_s
                do k=0,nz_s
                    do i=0,nxm
                        u2(i,k,j)=tmp(i+1,j+1,k+1)
                    end do
                end do
            end do
        case (3)
            do j=0,ny_s
                do k=0,nz_s
                    do i=0,nxm
                        u3(i,k,j)=tmp(i+1,j+1,k+1)
                    end do
                end do
            end do
        end select

! Close dataset
        call h5sclose_f(filspace_id, error)
        call h5dclose_f(dset_id, error)

    else if (.not.create_new_th(max(1,ith-3))) then
! Check to make sure that we should read in this scalar
        call h5dopen_f(file_id,trim(dname),dset_id,error)
        call h5dget_space_f(dset_id,filspace_id_th,error)
! Select hyperslab in the file.
        call h5sselect_hyperslab_f (filspace_id_th, H5S_SELECT_SET_F, &
                        offset_th, count, error, stride, block)
        call h5sselect_hyperslab_f (memspace_id_th, H5S_SELECT_SET_F, &
                        offset_m, count, error, stride, block)
! Read the dataset collectively
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmp_th, &
                dimsm_th, error, file_space_id = filspace_id_th, &
                mem_space_id = memspace_id_th, xfer_prp = plist_id_w)
        do j=0,ny_s_th
            do k=0,nz_s_th
                do i=0,nxm_th
                    th(i,k,j,ith-3)=tmp_th(i+1,j+1,k+1)
                end do
            end do
        end do

! Close dataset
        call h5sclose_f(filspace_id_th, error)
        call h5dclose_f(dset_id, error)
    end if

end do

! Decide whether to compute the pressure or to read
call h5lexists_f(file_id, 'p', read_pressure, error)
if (read_pressure) then
    dname="p"
    call h5dopen_f(file_id,trim(dname),dset_id,error)
    call h5dget_space_f(dset_id,filspace_id,error)

! Select hyperslab in the file.
    call h5sselect_hyperslab_f (filspace_id, h5s_select_set_f, &
                        offset_f, count, error, stride, block)
    call h5sselect_hyperslab_f (memspace_id, h5s_select_set_f, &
                        offset_m, count, error, stride, block)

! Write the dataset collectively
    call h5dread_f(dset_id, h5t_native_double, tmp, &
                dimsm, error, file_space_id = filspace_id, &
                mem_space_id = memspace_id, xfer_prp = plist_id_w)

    do j=0,ny_s
        do k=0,nz_s
            do i=0,nxm
                p(i,k,j)=tmp(i+1,j+1,k+1)
            end do
        end do
    end do
! Close dataset
    call h5sclose_f(filspace_id, error)
    call h5dclose_f(dset_id, error)
    call fft_xzy_mpi_to_fourier(p,cp)
end if

! Close the dataspace for the memory
call h5sclose_f(memspace_id, error)
call h5sclose_f(memspace_id_th, error)

! Close the properties for the reading
call h5pclose_f(plist_id_w, error)

! Close file
call h5fclose_f(file_id, error)
call h5close_f(error)

if (variable_dt) then
call courant_mpi
end if

! Convert to Fourier space
call fft_xzy_mpi_to_fourier(u1,cu1)
call fft_xzy_mpi_to_fourier(u2,cu2)
call fft_xzy_mpi_to_fourier(u3,cu3)
do ith=1,n_th
    if (.not.create_new_th(ith)) then
        s1(:,:,:)=th(:,:,:,ith)
        call fft_xzy_mpi_to_fourier_th(s1,cs1)
        cth(:,:,:,ith)=cs1(:,:,:)
    end if
end do

end subroutine ReadHDF5

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine init_movie
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use hdf5

include 'header'

character*55 fname

!     HDF5 ------------------------------------------------------

!     Identifiers
integer(hid_t) :: file_id, gid, plist_id_d, i

integer :: arank, error, comm_mpi, ranky_mov, rankz_mov
integer(hsize_t),dimension(2) :: adims
integer(size_t)               :: tdim
integer(hid_t)                :: aid,tspace,ttype
integer, dimension(1)         :: tint(6)
character*20                  :: sttimec


! Initialize interface
call h5open_f(error)

do i=1,3
    select case(i)
        case(1)
        fname='movie_xy.h5'
        comm_mpi=mpi_comm_y
    case(2)
        fname='movie_xz.h5'
        comm_mpi=mpi_comm_z
    case(3)
        fname='movie_yz.h5'
        comm_mpi=mpi_comm_world
    end select
    ranky_mov = ny_mov/(ny_s+1)
    rankz_mov = nz_mov/(nz_s+1)
    if ( ((i == 1).and.(rankz == rankz_mov)) .or. &
            ((i == 2).and.(ranky == ranky_mov)) .or. &
            (i == 3) ) then
! Setup file access property list with parallel I/O access
        call h5pcreate_f(h5p_file_access_f, plist_id_d, error)
        call h5pset_fapl_mpio_f(plist_id_d, comm_mpi, &
                    mpi_info_null, error)

! Create the file collectively
        call h5fcreate_f(fname, h5f_acc_trunc_f, &
                            file_id, error, access_prp = plist_id_d)
        call h5pclose_f(plist_id_d, error)

! Write file attributes:

        ! ----------------------------
        ! Resolution
        adims=(/3,2/)
        arank=2
        call h5screate_simple_f(arank,adims,tspace, error)
        call h5acreate_f(file_id,'Resolution',H5T_STD_I32LE,tspace, aid, error)
        tint(1)=nx
        tint(2)=ny
        tint(3)=nz
        tint(4)=nx_th
        tint(5)=ny_th
        tint(6)=nz_th
        call h5awrite_f(aid,h5t_native_integer,tint,adims,error)
        call h5aclose_f(aid, error)
        ! Close dataspace
        call h5sclose_f(tspace, error)

        ! -----------------------------
        ! Date
        adims=(/1,20/)
        arank=1
        tdim=20
        call h5tcopy_f(H5T_FORTRAN_S1, ttype, error)
        call h5tset_size_f(ttype,tdim,error)
        call h5screate_simple_f(arank,adims,tspace, error)

        call h5acreate_f(file_id,'Date',ttype,tspace,aid,
        &                 error)
        call time_string(sttimec)
        call h5awrite_f(aid,ttype,trim(sttimec),adims,
        &                  error)
        call h5aclose_f(aid, error)
        call h5sclose_f(tspace,error)
        call h5tclose_f(ttype,error)
        ! -----------------------------

! Create the groups for each plane (with attribute)
        call h5screate_f(h5s_scalar_f,tspace,error)
        call h5acreate_f(file_id,'Samples',h5t_std_i32le,tspace,aid,error)
        call h5awrite_f(aid,h5t_native_integer,0,adims,error)
        call h5aclose_f(aid,error)
        call h5sclose_f(tspace,error)

! Close file and interface
        call h5fclose_f(file_id,error)
    end if
end do
call h5close_f(error)

end subroutine INIT_MOVIE

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine WriteHDF5_xyplane(iz)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use hdf5

include 'header'

real*8 tmp(nx,ny_s+1)
real*8 tmp_th(nx_th,ny_s_th+1)
!
!     HDF5 ------------------------------------------------------
!
!     Dataset names
character(len=10) :: gname, dname

!     Identifiers
integer(hid_t) :: file_id, dset_id
integer(hid_t) :: filspace_id, memspace_id
integer(hid_t) :: filspace_id_th, memspace_id_th

!     Identifiers
integer(hid_t) :: gt_id, selspace_id
integer(hid_t) :: plist_id_w,plist_id_d,plist_id_d_th

!     Dimensions in the memory and in the file
integer(hsize_t), dimension(2) :: dimsm,dimsf
integer(hsize_t), dimension(2) :: dimsm_th,dimsf_th

integer(hsize_t), dimension(2) :: chunk_dims, count, offset_f
integer(hsize_t), dimension(2) :: stride, block, offset_m
integer(hsize_t), dimension(2) :: chunk_dims_th, offset_th
integer(hsize_t), dimension(2) :: block_th

integer :: rHDF5 = 2, arank, nsamp, IZ, i, j
logical flage

integer(hsize_t),dimension(2) :: adims
integer(size_t)               :: tdim
integer(hid_t)                :: aid,tspace,ttype
integer, dimension(1)         :: tint(6)
character*80                  :: namnbuf
character*20                  :: sttimec

integer error, ith

! #### DEFINE WRITING PARAMETERS ####
dimsm = (/nx,ny_s+1/)
dimsf = (/nx,ny/)
block = dimsm
dimsm_th = (/nx_th,ny_s_th+1/)
dimsf_th = (/nx_th,ny_th/)
block_th = dimsm_th

chunk_dims = (/nx,1/)
chunk_dims_th = (/nx_th,1/)

! Stride and count for number of rows and columns in each dimension
stride = 1
count  = 1

! Offset determined by the rank of a processor
offset_f = (/ 0, ranky*(ny_s+1) /)
offset_th = (/ 0, ranky*(ny_s_th+1) /)
offset_m = 0

! #### SET UP GROUP STRUCTURE ####

! Initialize interface
call h5open_f(error)

! Setup file access property list with parallel I/O access
call h5pcreate_f(h5p_file_access_f, plist_id_d, error)
call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_y, mpi_info_null, error)

! Open the file collectively and the relevant group
call h5fopen_f('movie_xy.h5', H5F_ACC_RDWR_F, &
                file_id, error, access_prp = plist_id_d)
call h5pclose_f(plist_id_d, error)

! Check and update the number of samples
adims=1
call h5aopen_f(file_id,'Samples',aid,error)
call h5aread_f(aid,h5t_native_integer,nsamp,adims,error)
! Create group for this timestep
write(gname,'(1i0.4)') nsamp
call h5gcreate_f(file_id,gname,gt_id,error)
nsamp=nsamp+1
call h5awrite_f(aid,h5t_native_integer,nsamp,adims,error)
call h5aclose_f(aid,error)

! Write Time, Timestep, missing coordinate attributes
call h5screate_f(H5S_SCALAR_F,tspace,error)

call h5acreate_f(gt_id,'Time',H5T_IEEE_F64LE,tspace,
&                 aid, error)
call h5awrite_f(aid,H5T_NATIVE_DOUBLE,time,adims,error)
call h5aclose_f(aid, error)

call h5acreate_f(gt_id,'Timestep',H5T_STD_I32LE,tspace,
&                 aid, error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,time_step,adims,error)
call h5aclose_f(aid, error)

call h5acreate_f(gt_id,'z',H5T_IEEE_F64LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_DOUBLE,gz(iz),adims,error)
call h5aclose_f(aid,error)

call h5sclose_f(tspace,error)
! ----------------------------

! #### WRITE DATASETS ####

! Create property list for the chunked dataset creation
call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id_d, error)
call h5pset_chunk_f(plist_id_d, rHDF5, chunk_dims, error)
call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id_d_th, error)
call h5pset_chunk_f(plist_id_d_th, rHDF5, chunk_dims_th, error)

! Create the dataspaces for velocities and scalars
call h5screate_simple_f(rHDF5, dimsf, filspace_id, error)
call h5screate_simple_f(rHDF5, dimsm, memspace_id, error)
call h5screate_simple_f(rHDF5, dimsf_th, filspace_id_th, error)
call h5screate_simple_f(rHDF5, dimsm_th, memspace_id_th, error)

! Create property list for collective dataset write
call h5pcreate_f(h5p_dataset_xfer_f, plist_id_w, error)
call h5pset_dxpl_mpio_f(plist_id_w, h5fd_mpio_collective_f, error)
if (rankz==nz_mov/(nz_s+1)) then
    do ith=1,3+n_th

        select case(ith)
        case (1)
            do i=0,nxm
                do j=0,ny_s
                    tmp(i+1,j+1) = u1(i,iz,j)
                end do
            end do
            dname="u1"
        case (2)
            do i=0,nxm
                do j=0,ny_s
                    tmp(i+1,j+1) = u2(i,iz,j)
                end do
            end do
            dname="u2"
        case (3)
            do i=0,nxm
                do j=0,ny_s
                    tmp(i+1,j+1) = u3(i,iz,j)
                end do
            end do
            dname="u3"
        case (4:)
            do i=0,nxm_th
                do j=0,ny_s_th
                    tmp_th(i+1,j+1) = th(i,iz,j,ith-3)
                end do
            end do
            dname="th"//char(ith+45)
        end select
!      if (ith == 1) write(*,*) 'RANK ',RANK,U1(:,IZ,NY_S)
        if (ith <= 3) then
! If saving velocity, use velocity grid
            call h5dcreate_f(gt_id, trim(dname), H5T_IEEE_F64LE, &
                    filspace_id, dset_id, error, dcpl_id = plist_id_d)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
            call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F, &
                offset_f, count, error, stride, block)

            call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, &
                offset_m, count, error, stride, block)

! Write the dataset collectively
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                            tmp, &
                            dimsm, error, file_space_id = filspace_id, &
                            mem_space_id = memspace_id, xfer_prp = plist_id_w)

! Close dateset
            call h5dclose_f(dset_id, error)
        else
! If saving scalar, use TH grid
            call h5dcreate_f(gt_id, trim(dname), H5T_IEEE_F64LE, &
                    filspace_id_th, dset_id, error, dcpl_id = plist_id_d_th)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
            call h5sselect_hyperslab_f (filspace_id_th, h5s_select_set_f, &
                        offset_th, count, error, stride, block_th)

            call h5sselect_hyperslab_f (memspace_id_th, h5s_select_set_f, &
                        offset_m, count, error, stride, block_th)

! Write the dataset collectively
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                            tmp_th, &
                            dimsm_th, error, file_space_id = filspace_id_th, &
                            mem_space_id = memspace_id_th, xfer_prp = plist_id_w)

! Close dataset
            call h5dclose_f(dset_id, error)
        end if
    end do
end if

!     Close the dataspace for the memory and for the file
call h5sclose_f(filspace_id, error)
call h5sclose_f(memspace_id, error)
call h5sclose_f(filspace_id_th, error)
call h5sclose_f(memspace_id_th, error)

!     Close the properties for the dataspace creation and the writing
call h5pclose_f(plist_id_d, error)
call h5pclose_f(plist_id_d_th, error)
call h5pclose_f(plist_id_w, error)

! Close groups
call h5gclose_f(gt_id,error)
call h5fclose_f(file_id, error)
call h5close_f(error)

! call mpi_finalize(ierror)
! stop

end subroutine WriteHDF5_xyplane

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine WriteHDF5_xzplane(iy)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use hdf5

include 'header'


real*8 tmp(nx,nz_s+1)
real*8 tmp_th(nx_th,nz_s_th+1)
!
!     HDF5 ------------------------------------------------------
!
!     Dataset names
character(len=10) :: gname, dname

!     Identifiers
integer(hid_t) :: file_id, dset_id
integer(hid_t) :: filspace_id, memspace_id
integer(hid_t) :: filspace_id_th, memspace_id_th

!     Identifiers
integer(hid_t) :: gt_id, selspace_id
integer(hid_t) :: plist_id_w,plist_id_d,plist_id_d_th

!     Dimensions in the memory and in the file
integer(hsize_t), dimension(2) :: dimsm,dimsf
integer(hsize_t), dimension(2) :: dimsm_th,dimsf_th

integer(hsize_t), dimension(2) :: chunk_dims, count, offset_f
integer(hsize_t), dimension(2) :: stride, block, offset_m
integer(hsize_t), dimension(2) :: chunk_dims_th, offset_th
integer(hsize_t), dimension(2) :: block_th

integer :: rHDF5 = 2, arank = 1, nsamp, iy, i, k
logical flage

integer(hsize_t),dimension(1) :: adims
integer(size_t)               :: tdim
integer(hid_t)                :: aid,tspace,ttype
integer, dimension(1)         :: tint(6)
character*80                  :: namnbuf
character*20                  :: sttimec

integer error, ith

! #### DEFINE WRITING PARAMETERS ####
dimsm = (/nx,nz_s+1/)
dimsf = (/nx,nz/)
block = dimsm
dimsm_th = (/nx_th,nz_s_th+1/)
dimsf_th = (/nx_th,nz_th/)
block_th = dimsm_th

chunk_dims = (/nx,nz_s+1/)
chunk_dims_th = (/nx_th,nz_s_th+1/)

! Stride and count for number of rows and columns in each dimension
stride = 1
count  = 1

! Offset determined by the rank of a processor
offset_f = (/ 0, rankz*(nz_s+1) /)
offset_th = (/ 0, rankz*(nz_s_th+1) /)
offset_m = 0

! #### SET UP GROUP STRUCTURE ####

! Initialize interface
call h5open_f(error)

! Setup file access property list with parallel I/O access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id_d, error)
call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_z, mpi_info_null, error)

! Open the file collectively and the relevant group
call h5fopen_f('movie_xz.h5', H5F_ACC_RDWR_F, &
                    file_id, error, access_prp = plist_id_d)
call h5pclose_f(plist_id_d, error)

! Check and update the number of samples
call h5aopen_f(file_id,'Samples',aid,error)
call h5aread_f(aid,H5T_NATIVE_INTEGER,nsamp,adims,error)
! Create group for this timestep
write(gname,'(1i0.4)') nsamp
call h5gcreate_f(file_id,gname,gt_id,error)
nsamp=nsamp+1
call h5awrite_f(aid,h5t_native_integer,nsamp,adims,error)
call h5aclose_f(aid,error)

! Write Time, Timestep, missing coordinate attributes
call h5screate_f(H5S_SCALAR_F,tspace,error)

call h5acreate_f(gt_id,'Time',H5T_IEEE_F64LE,tspace, aid, error)
call h5awrite_f(aid,H5T_NATIVE_DOUBLE,time,adims,error)
call h5aclose_f(aid, error)

call h5acreate_f(gt_id,'Timestep',H5T_STD_I32LE,tspace, aid, error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,time_step,adims,error)
call h5aclose_f(aid, error)

call h5acreate_f(gt_id,'y',h5t_ieee_f64le,tspace,aid,error)
call h5awrite_f(aid,h5t_native_double,gy(iy),adims,error)
call h5aclose_f(aid,error)

call h5sclose_f(tspace,error)
! ----------------------------

! #### WRITE DATASETS ####

! Create property list for the chunked dataset creation
call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id_d, error)
call h5pset_chunk_f(plist_id_d, rHDF5, chunk_dims, error)
call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id_d_th, error)
call h5pset_chunk_f(plist_id_d_th, rHDF5, chunk_dims_th, error)

! Create the dataspaces for velocities and scalars
call h5screate_simple_f(rHDF5, dimsf, filspace_id, error)
call h5screate_simple_f(rHDF5, dimsm, memspace_id, error)
call h5screate_simple_f(rHDF5, dimsf_th, filspace_id_th, error)
call h5screate_simple_f(rHDF5, dimsm_th, memspace_id_th, error)

! Create property list for collective dataset write
call h5pcreate_f(h5p_dataset_xfer_f, plist_id_w, error)
call h5pset_dxpl_mpio_f(plist_id_w, h5fd_mpio_collective_f, error)
if (ranky==ny_mov/(ny_s+1)) then
    do ith=1,3+n_th

        select case(ith)
        case (1)
            do i=0,nxm
                do k=0,nz_s
                    tmp(i+1,k+1) = u1(i,k,iy)
                end do
            end do
            dname="u1"
        case (2)
            do i=0,nxm
                do k=0,nz_s
                    tmp(i+1,k+1) = u2(i,k,iy)
                end do
            end do
            dname="u2"
        case (3)
            do i=0,nxm
                do k=0,nz_s
                    tmp(i+1,k+1) = u3(i,k,iy)
                end do
            end do
            dname="u3"
        case (4:)
            do i=0,nxm_th
                do k=0,nz_s_th
                    tmp_th(i+1,k+1) = th(i,k,iy,ith-3)
                end do
            end do
            dname="th"//char(ith+45)
        end select
        if (ith <= 3) then
! If saving velocity, use velocity grid
            call h5dcreate_f(gt_id, trim(dname), H5T_IEEE_F64LE, &
                    filspace_id, dset_id, error, dcpl_id = plist_id_d)

! Select hyperslab in the file.
    !     call h5dget_space_f(dsetur_id, selspace_id, error)
            call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F, &
                        offset_f, count, error, stride, block)

            call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, &
                        offset_m, count, error, stride, block)

! Write the dataset collectively
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                    tmp, &
                    dimsm, error, file_space_id = filspace_id, &
                    mem_space_id = memspace_id, xfer_prp = plist_id_w)

! Close dateset
            call h5dclose_f(dset_id, error)
        else
! If saving scalar, use TH grid
            call h5dcreate_f(gt_id, trim(dname), H5T_IEEE_F64LE, &
                    filspace_id_th, dset_id, error, dcpl_id = plist_id_d_th)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
            call h5sselect_hyperslab_f (filspace_id_th, H5S_SELECT_SET_F, &
                    offset_th, count, error, stride, block_th)

            call h5sselect_hyperslab_f (memspace_id_th, H5S_SELECT_SET_F, &
                    offset_m, count, error, stride, block_th)

! Write the dataset collectively
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                    tmp_th, &
                    dimsm_th, error, file_space_id = filspace_id_th, &
                    mem_space_id = memspace_id_th, xfer_prp = plist_id_w)

! Close dataset
            call h5dclose_f(dset_id, error)
        end if
    end do
end if

!     Close the dataspace for the memory and for the file
call h5sclose_f(filspace_id, error)
call h5sclose_f(memspace_id, error)
call h5sclose_f(filspace_id_th, error)
call h5sclose_f(memspace_id_th, error)

!     Close the properties for the dataspace creation and the writing
call h5pclose_f(plist_id_d, error)
call h5pclose_f(plist_id_d_th, error)
call h5pclose_f(plist_id_w, error)

! Close groups
call h5gclose_f(gt_id,error)
call h5fclose_f(file_id, error)
call h5close_f(error)

! call mpi_finalize(ierror)
! stop

end subroutine WriteHDF5_xzplane



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine WriteHDF5_yzplane(ix)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use hdf5

include 'header'

real*8 tmp(ny_s+1,nz_s+1)
real*8 tmp_th(ny_s_th+1,nz_s_th+1)
!
!     HDF5 ------------------------------------------------------
!
!     Dataset names
character(len=10) :: gname, dname

!     Identifiers
integer(hid_t) :: file_id, dset_id
integer(hid_t) :: filspace_id, memspace_id
integer(hid_t) :: filspace_id_th, memspace_id_th

!     Identifiers
integer(hid_t) :: gt_id, selspace_id
integer(hid_t) :: plist_id_w,plist_id_d,plist_id_d_th

!     Dimensions in the memory and in the file
integer(hsize_t), dimension(2) :: dimsm,dimsf
integer(hsize_t), dimension(2) :: dimsm_th,dimsf_th

integer(hsize_t), dimension(2) :: chunk_dims, count, offset_f
integer(hsize_t), dimension(2) :: stride, block, offset_m
integer(hsize_t), dimension(2) :: chunk_dims_th, offset_th
integer(hsize_t), dimension(2) :: block_th

integer :: rHDF5 = 2, arank = 1, nsamp, IX, j, k
logical flage

integer(hsize_t),dimension(1) :: adims
integer(size_t)               :: tdim
integer(hid_t)                :: aid,tspace,ttype
integer, dimension(1)         :: tint(6)
character*80                  :: namnbuf
character*20                  :: sttimec

integer error, ith

! #### DEFINE WRITING PARAMETERS ####
dimsm = (/ny_s+1,nz_s+1/)
dimsf = (/ny,nz/)
block = dimsm
dimsm_th = (/ny_s_th+1,nz_s_th+1/)
dimsf_th = (/ny_th,nz_th/)
block_th = dimsm_th

chunk_dims = (/1,nz_s+1/)
chunk_dims_th = (/1,nz_s_th+1/)

! Stride and count for number of rows and columns in each dimension
stride = 1
count  = 1

! Offset determined by the rank of a processor
offset_f = (/ ranky*(ny_s+1), rankz*(nz_s+1) /)
offset_th = (/ ranky*(ny_s_th+1), rankz*(nz_s_th+1) /)
offset_m = 0
! #### SET UP GROUP STRUCTURE ####

! Initialize interface
call h5open_f(error)

! Setup file access property list with parallel I/O access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id_d, error)
call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_world, mpi_info_null, error)
! Open the file collectively and the relevant group
call h5fopen_f('movie_yz.h5', H5F_ACC_RDWR_F, &
                    file_id, error, access_prp = plist_id_d)
call h5pclose_f(plist_id_d, error)

! Check and update the number of samples
call h5aopen_f(file_id,'Samples',aid,error)
call h5aread_f(aid,H5T_NATIVE_INTEGER,nsamp,adims,error)
! Create group for this timestep
write(gname,'(1i0.4)') nsamp
call h5gcreate_f(file_id,gname,gt_id,error)
nsamp=nsamp+1
call h5awrite_f(aid,H5T_NATIVE_INTEGER,nsamp,adims,error)
call h5aclose_f(aid,error)
! Write Time, Timestep, missing coordinate attributes
call h5screate_f(H5S_SCALAR_F,tspace,error)

call h5acreate_f(gt_id,'Time',H5T_IEEE_F64LE,tspace, aid, error)
call h5awrite_f(aid,H5T_NATIVE_DOUBLE,TIME,adims,error)
call h5aclose_f(aid, error)

call h5acreate_f(gt_id,'Timestep',H5T_STD_I32LE,tspace, aid, error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,TIME_STEP,adims,error)
call h5aclose_f(aid, error)
call h5acreate_f(gt_id,'x',H5T_IEEE_F64LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_DOUBLE,GX(IX),adims,error)
call h5aclose_f(aid,error)

call h5sclose_f(tspace,error)
! ----------------------------
call mpi_barrier(mpi_comm_world,ierror)
! #### WRITE DATASETS ####

! Create property list for the chunked dataset creation
call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id_d, error)
call h5pset_chunk_f(plist_id_d, rHDF5, chunk_dims, error)
call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id_d_th, error)
call h5pset_chunk_f(plist_id_d_th, rHDF5, chunk_dims_th, error)

! Create the dataspaces for velocities and scalars
call h5screate_simple_f(rHDF5, dimsf, filspace_id, error)
call h5screate_simple_f(rHDF5, dimsm, memspace_id, error)
call h5screate_simple_f(rHDF5, dimsf_th, filspace_id_th, error)
call h5screate_simple_f(rHDF5, dimsm_th, memspace_id_th, error)

! Create property list for collective dataset write
call h5pcreate_f(h5p_dataset_xfer_f, plist_id_w, error)
call h5pset_dxpl_mpio_f(plist_id_w, h5fd_mpio_collective_f, error)

do ith=1,3+n_th

    select case(ith)
    case (1)
        do j=0,ny_s
            do k=0,nz_s
                tmp(j+1,k+1) = u1(ix,k,j)
            end do
        end do
        dname="u1"
    case (2)
        do j=0,ny_s
            do k=0,nz_s
                tmp(j+1,k+1) = u2(ix,k,j)
            end do
        end do
        dname="u2"
    case (3)
        do j=0,ny_s
            do k=0,nz_s
                tmp(j+1,k+1) = u3(ix,k,j)
            end do
        end do
        dname="u3"
    case (4:)
        do j=0,ny_s_th
            do k=0,nz_s_th
                tmp_th(j+1,k+1) = th(ix,k,j,ith-3)
            end do
        end do
        dname="th"//char(ith+45)
    end select
    if (ith <= 3) then
! If saving velocity, use velocity grid
        call h5dcreate_f(gt_id, trim(dname), H5T_IEEE_F64LE, &
                    filspace_id, dset_id, error, dcpl_id = plist_id_d)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
        call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F, &
                    offset_f, count, error, stride, block)

        call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, &
                    offset_m, count, error, stride, block)

! Write the dataset collectively
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                        tmp, &
                        dimsm, error, file_space_id = filspace_id, &
                        mem_space_id = memspace_id, xfer_prp = plist_id_w)

! Close dateset
        call h5dclose_f(dset_id, error)
    else
! If saving scalar, use TH grid
        call h5dcreate_f(gt_id, trim(dname), H5T_IEEE_F64LE, &
                filspace_id_th, dset_id, error, dcpl_id = plist_id_d_th)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
        call h5sselect_hyperslab_f (filspace_id_th, H5S_SELECT_SET_F, &
                    offset_th, count, error, stride, block_th)

        call h5sselect_hyperslab_f (memspace_id_th, H5S_SELECT_SET_F, &
                    offset_m, count, error, stride, block_th)

! Write the dataset collectively
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                        tmp_th, &
                        dimsm_th, error, file_space_id = filspace_id_th, &
                        mem_space_id = memspace_id_th, xfer_prp = plist_id_w)

! Close dataset
        call h5dclose_f(dset_id, error)
    end if
end do

!     Close the dataspace for the memory and for the file
call h5sclose_f(filspace_id, error)
call h5sclose_f(memspace_id, error)
call h5sclose_f(filspace_id_th, error)
call h5sclose_f(memspace_id_th, error)

!     Close the properties for the dataspace creation and the writing
call h5pclose_f(plist_id_d, error)
call h5pclose_f(plist_id_d_th, error)
call h5pclose_f(plist_id_w, error)

      ! Close groups
call h5gclose_f(gt_id,error)
call h5fclose_f(file_id, error)
call h5close_f(error)

! call mpi_finalize(ierror)
! stop

end subroutine WriteHDF5_yzplane


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine init_stats
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use hdf5

include 'header'

!     HDF5 ------------------------------------------------------

!     Identifiers
integer(hid_t) :: file_id, gid, plist_id_d

integer :: arank, error
integer(hsize_t),dimension(2) :: adims
integer(size_t)               :: tdim
integer(hid_t)                :: aid,tspace,ttype
character*20                  :: sttimec


! Initialize interface
call h5open_f(error)

! Setup file access property list with parallel I/O access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id_d, error)
call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_world, &
                mpi_info_null, error)

! Create the file collectively
call h5fcreate_f('stats.h5', H5F_ACC_TRUNC_F, &
                    file_id, error, access_prp = plist_id_d)
call h5pclose_f(plist_id_d, error)

! Write file attributes:

! -----------------------------
! Date
adims=(/1,20/)
arank=1
tdim=20
call h5tcopy_f(H5T_FORTRAN_S1, ttype, error)
call h5tset_size_f(ttype,tdim,error)
call h5screate_simple_f(arank,adims,tspace, error)

call h5acreate_f(file_id,'Date',ttype,tspace,aid, error)
call time_string(sttimec)
call h5awrite_f(aid,ttype,trim(sttimec),adims, error)
call h5aclose_f(aid, error)
call h5sclose_f(tspace,error)
call h5tclose_f(ttype,error)
! -----------------------------

! Create the groups for each quantity (with attribute)
call h5gcreate_f(file_id,'U1rms',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'U2rms',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'U3rms',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'THrms',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'THflux',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'epsilon',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'chi',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

! Close file and interface
call h5fclose_f(file_id,error)
call h5close_f(error)

end subroutine init_stats


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine WriteStatH5(gname,stat)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use hdf5

include 'header'

real*8 stat
!
!     HDF5 ------------------------------------------------------
!
!     Dataset names
character(len=10) :: gname, dname

!     Identifiers
integer(hid_t) :: file_id, dset_id
integer(hid_t) :: filspace_id, memspace_id
integer(hid_t) :: filspace_id_th, memspace_id_th

!     Identifiers
integer(hid_t) :: gid, gt_id, selspace_id
!      integer(hid_t) :: plist_id_w,plist_id_d,plist_id_d_th

!     Dimensions in the memory and in the file
integer(hsize_t), dimension(1) :: dims

integer :: rHDF5 = 1, arank = 1, nsamp
logical flage

integer(hsize_t),dimension(1) :: adims = 1
integer(size_t)               :: tdim
integer(hid_t)                :: aid,tspace,ttype
character*80                  :: namnbuf
character*20                  :: sttimec

integer error, ith

! #### DEFINE WRITING PARAMETERS ####
dims = 1

! #### SET UP GROUP STRUCTURE ####

! Initialize interface
call h5open_f(error)

! Open the file collectively and the relevant group
call h5fopen_f('stats.h5', H5F_ACC_RDWR_F,file_id, error)
call h5gopen_f(file_id,'/'//trim(gname),gid,error)

! Check and update the number of samples
call h5aopen_f(gid,'Samples',aid,error)
call h5aread_f(aid,H5T_NATIVE_INTEGER,nsamp,adims,error)
! Create dataset name for this timestep
write(dname,'(1i0.4)') nsamp
nsamp=nsamp+1
call h5awrite_f(aid,H5T_NATIVE_INTEGER,nsamp,adims,error)
call h5aclose_f(aid,error)

call h5screate_f(H5S_SCALAR_F,tspace,error)

! Write statistic to file
call h5dcreate_f(gid,dname,H5T_IEEE_F64LE,tspace,dset_id,error)
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,Stat,dims,error)

! Write Time, Timestep attributes
call h5acreate_f(dset_id,'Time',H5T_IEEE_F64LE,tspace, aid, error)
call h5awrite_f(aid,H5T_NATIVE_DOUBLE,TIME,adims,error)
call h5aclose_f(aid, error)

call h5acreate_f(dset_id,'Timestep',H5T_STD_I32LE,tspace, aid, error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,TIME_STEP,adims,error)
call h5aclose_f(aid, error)


call h5dclose_f(dset_id,error)

call h5sclose_f(tspace,error)
! ----------------------------

! Close groups
call h5gclose_f(gid,error)
call h5fclose_f(file_id, error)
call h5close_f(error)

end subroutine WriteStatH5


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine init_mean
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use hdf5

include 'header'

!     Identifiers
integer(hid_t) :: file_id, gid, plist_id_d

integer :: arank, error
integer(hsize_t), dimension(1) :: adims
integer(hid_t)                 :: aid, tspace

! Initialize interface
call h5open_f(error)

! Setup file access property list with parallel I/O access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id_d, error)
call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_y, &
                mpi_info_null, error)

! Create the file collectively
call h5fcreate_f('mean.h5', H5F_ACC_TRUNC_F, &
                    file_id, error, access_prp = plist_id_d)
call h5pclose_f(plist_id_d, error)

! Create the groups for each quantity (with attributes)
call h5gcreate_f(file_id,'U1me',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'U3me',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'THme',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'epsilon',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'chi',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'U1U2',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'U3U2',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'THflux',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'U1rms',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'U2rms',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'U3rms',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'THrms',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

! Close file and interface
call h5fclose_f(file_id,error)
call h5close_f(error)

end subroutine

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine WriteMeanH5(gname,Ume)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use HDF5

include 'header'

real*8 Ume(0:ny_s)

!     Dataset names
character(len=10) :: gname, dname

!     Identifiers
integer(hid_t) :: file_id, gid, aid, dset_id
integer(hid_t) :: filspace_id, memspace_id
integer(hid_t) :: plist_id_d, plist_id_w

!     Dimensions
integer(hsize_t), dimension(1) :: adims, dimsm, dimsf
integer(hsize_t), dimension(1) :: offset_f, offset_m
integer(hsize_t), dimension(1) :: stride, block, count
integer :: rHDF5 = 1

integer error, nsamp

! #### DEFINE WRITING PARAMETERS ####
dimsm=ny_s+1
dimsf=ny
block=dimsm
stride=1
count=1      
offset_f = ranky*(ny_s+1)
offset_m = 0

! #### SET UP GROUP STRUCTURE ####

! Initialize interface
call h5open_f(error)

! Setup file access property list with parallel I/O access
call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id_d, error)
call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_y, mpi_info_null, error)

! Open the file collectively and the relevant group
call h5fopen_f('mean.h5', H5F_ACC_RDWR_F, &
                    file_id, error, access_prp = plist_id_d)
call h5pclose_f(plist_id_d, error)
call h5gopen_f(file_id,'/'//trim(gname),gid,error)

! Check number of samples
adims=1
call h5aopen_f(gid,'Samples',aid,error)
call h5aread_f(aid,H5T_NATIVE_INTEGER,nsamp,adims,error)
! Create dataset name for this timestep
write(dname,'(1i0.4)') nsamp
nsamp=nsamp+1
! Update number of samples
call h5awrite_f(aid,H5T_NATIVE_INTEGER,nsamp,adims,error)
call h5aclose_f(aid,error)

! #### WRITE DATASETS ####

! Create property list for the dataset creation
call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id_d, error)

! Create the dataspace
call h5screate_simple_f(rHDF5, dimsf, filspace_id, error)
call h5screate_simple_f(rHDF5, dimsm, memspace_id, error)

! Create property list for collective dataset write
call h5pcreate_f(H5P_DATASET_XFER_F, plist_id_w, error)
call h5pset_dxpl_mpio_f(plist_id_w, H5FD_MPIO_COLLECTIVE_F, error)

! Create dataset
call h5dcreate_f(gid,dname,H5T_IEEE_F64LE,filspace_id, &
                    dset_id,error, dcpl_id = plist_id_d)

! Select hyperslab in the file.
call h5sselect_hyperslab_f(filspace_id, H5S_SELECT_SET_F, &
                offset_f, count, error, stride, block)
call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F, &
                offset_m, count, error, stride, block)

! Write the dataset collectively
call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                Ume, &
                dimsm, error, file_space_id = filspace_id, &
                mem_space_id = memspace_id, xfer_prp = plist_id_w)

! Close dataset, spaces, properties, group, file & interface
call h5dclose_f(dset_id,error)
call h5sclose_f(filspace_id,error)
call h5sclose_f(memspace_id,error)
call h5pclose_f(plist_id_d,error)
call h5pclose_f(plist_id_w,error)
call h5gclose_f(gid,error)
call h5fclose_f(file_id,error)
call h5close_f(error)

end subroutine


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine init_spectra
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use hdf5

include 'header'

!     HDF5 ------------------------------------------------------

!     Identifiers
integer(hid_t) :: file_id, gid

integer :: arank, error
integer(hsize_t), dimension(1)  :: adims
integer(hid_t)    :: aid,tspace
character*20      :: sttimec


! Initialize interface
call h5open_f(error)

! Create the file collectively
call h5fcreate_f('spectra.h5', H5F_ACC_TRUNC_F, file_id, error)

! Create the groups for each quantity (with attribute)
call h5gcreate_f(file_id,'U1',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'U2',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'U3',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

call h5gcreate_f(file_id,'TH1',gid,error)
call h5screate_f(H5S_SCALAR_F,tspace,error)
call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
call h5aclose_f(aid,error)
call h5sclose_f(tspace,error)
call h5gclose_f(gid,error)

! Close file and interface
call h5fclose_f(file_id,error)
call h5close_f(error)

end subroutine init_spectra


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine WriteSpectrumH5(gname,spectrum)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
use hdf5

include 'header'

real*8 spectrum(0:tnky)
!
!     HDF5 ------------------------------------------------------
!
!     Dataset names
character(len=10) :: gname, dname

!     Identifiers
integer(hid_t) :: file_id, dset_id
integer(hid_t) :: gid

!     Dimensions in the memory and in the file
integer(hsize_t), dimension(1) :: dims

integer :: nsamp

integer(hsize_t),dimension(1) :: adims = 1
integer(hid_t)                :: aid,tspace
integer :: rHDF5 = 1
integer error

! #### DEFINE WRITING PARAMETERS ####
dims = tnky+1

! #### SET UP GROUP STRUCTURE ####

! Initialize interface
call h5open_f(error)

! Open the file collectively and the relevant group
call h5fopen_f('spectra.h5', H5F_ACC_RDWR_F,file_id, error)
call h5gopen_f(file_id,'/'//trim(gname),gid,error)

! Check and update the number of samples
call h5aopen_f(gid,'Samples',aid,error)
call h5aread_f(aid,H5T_NATIVE_INTEGER,nsamp,adims,error)
! Create dataset name for this timestep
write(dname,'(1i0.4)') nsamp
nsamp=nsamp+1
call h5awrite_f(aid,H5T_NATIVE_INTEGER,nsamp,adims,error)
call h5aclose_f(aid,error)

call h5screate_simple_f(rHDF5,dims,tspace,error)

! Write spectrum to file
call h5dcreate_f(gid,dname,H5T_IEEE_F64LE,tspace,dset_id,error)
call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,spectrum,dims,error)

call h5dclose_f(dset_id,error)

call h5sclose_f(tspace,error)
! ----------------------------

! Close groups
call h5gclose_f(gid,error)
call h5fclose_f(file_id, error)
call h5close_f(error)

end subroutine WriteSpectrumH5