C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ReadHDF5(FNAME)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      use hdf5

      INCLUDE 'header'

      CHARACTER*85 FNAME
      LOGICAL FINAL, READ_PRESSURE

      REAL*8 tmp(NX,NY,NZ)
      REAL*8 tmp_th(NX_TH,NY_TH,NZ_TH)

c     Dataset names
      character(len=10) :: dname

c     Identifiers
      integer(hid_t) :: file_id, dset_id
      integer(hid_t) :: filspace_id, memspace_id
      integer(hid_t) :: filspace_id_th, memspace_id_th

!     Identifiers
      integer(hid_t) :: selspace_id
      integer(hid_t) :: plist_id_w,plist_id_d

c     Dimensions in the memory and in the file
      integer(hsize_t), dimension(3) :: dimsm
      integer(hsize_t), dimension(3) :: dimsm_th

      integer(hsize_t), dimension(3) :: count, offset_f, offset_m
      integer(hsize_t), dimension(3) :: stride, block
      integer(hsize_t), dimension(3) :: offset_th, block_th

      integer :: rHDF5 = 3, i,j,k

      integer(hsize_t),dimension(2) :: adims
      integer(hid_t)                :: aid
      integer, dimension(2)         :: tint(3,2)
      character*80                  :: namnbuf
      character*20                  :: sttimec

      integer error, ith

      double precision En(4)

      dimsm(1:3) = (/NX,NY,NZ/)
      block = dimsm
      dimsm_th(1:3) = (/NX_TH,NY_TH,NZ_TH/)
      block_th = dimsm_th

!     Stride and count for number of rows and columns in each dimension
      stride = 1
      count  = 1

!     Offset determined by the rank of a processor
      offset_f(1:3) = 0
      offset_th(1:3) = 0
      offset_m(1:3) = 0

!     Initialize interface
      call h5open_f(error)

!     Open the file
      call h5fopen_f(trim(FNAME), H5F_ACC_RDONLY_F,
     &                 file_id, error)

!     -----------------------------
!     Resolution
!     -----------------------------
      adims=(/3,2/)
      call h5aopen_by_name_f(file_id,'.','Resolution',aid,error)
      call h5aread_f(aid,H5T_NATIVE_INTEGER,tint,adims,error)
      call h5aclose_f(aid, error)
      ! Check that the resolution is of the same kind
      if ((tint(1,1).ne. NX)  .or.  (tint(2,1).ne.NY)  .or.
     &    (tint(3,1).ne.NZ) .or. (tint(1,2).ne.NX_TH) .or.
     &    (tint(2,2).ne.NY_TH) .or. (tint(3,2).ne.NZ_TH)) then
        write(*,*) ' Error. File and program have ',
     &        'different resolutions. '
        write(*,*) ' Program: ', NX,NY,NZ,NX_TH,NY_TH,NZ_TH
        write(*,*) ' File   : ', tint(1:3,1),tint(1:3,2)
        stop
      end if

!     -----------------------------
!     Time stamps
!     -----------------------------
      adims=(/1,80/)
      call h5aopen_by_name_f(file_id,'.','Time',aid,error)
      call h5aread_f(aid,H5T_NATIVE_DOUBLE,TIME,adims,error)
      call h5aclose_f(aid, error)

      call h5aopen_by_name_f(file_id,'.','Timestep',aid,error)
      call h5aread_f(aid,H5T_NATIVE_INTEGER,TIME_STEP,adims,error)
      call h5aclose_f(aid, error)
!     -----------------------------

!     Dataspaces in memory
      call h5screate_simple_f(rHDF5, dimsm, memspace_id, error)
      call h5screate_simple_f(rHDF5, dimsm_th, memspace_id_th, error)

      do ith=1,3+N_TH
        write(*,*) 'ith: ',ith
!     Here it starts the loop--->
        select case(ith)
        case (1)
          dname="U"
        case (2)
          dname="V"
        case (3)
          dname="W"
        case (4:)
          dname="TH"//CHAR(ith+45)
        end select

        if (ith.le.3) then
          call h5dopen_f(file_id,trim(dname),dset_id,error)
          call h5dget_space_f(dset_id,filspace_id,error)
! Select hyperslab in the file.
          call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
     &        offset_f, count, error, stride, block)
          call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block)
! Read the dataset collectively
          call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmp,
     +        dimsm, error, file_space_id = filspace_id,
     +        mem_space_id = memspace_id)
          select case(ith)
          case (1)
          do j=0,NY-1
            do k=0,NZ-1
              do i=0,NX-1
                U1(i,k,j)=tmp(i+1,j+1,k+1)
                end do
            end do
          end do
          case (2)
          do j=0,NY-1
            do k=0,NZ-1
              do i=0,NX-1
                U2(i,k,j)=tmp(i+1,j+1,k+1)
                end do
            end do
          end do
          case (3)
          do j=0,NY-1
            do k=0,NZ-1
              do i=0,NX-1
                U3(i,k,j)=tmp(i+1,j+1,k+1)
                end do
            end do
          end do
          end select

! Close dataset
          call h5sclose_f(filspace_id, error)
          call h5dclose_f(dset_id, error)

        else
! Check to make sure that we should read in this scalar
          call h5dopen_f(file_id,trim(dname),dset_id,error)
          call h5dget_space_f(dset_id,filspace_id_th,error)
! Select hyperslab in the file.
          call h5sselect_hyperslab_f (filspace_id_th, H5S_SELECT_SET_F,
     &        offset_th, count, error, stride, block)
          call h5sselect_hyperslab_f (memspace_id_th, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block)
! Read the dataset collectively
          call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmp_th,
     +        dimsm_th, error, file_space_id = filspace_id_th,
     +        mem_space_id = memspace_id_th)
          do j=0,NY_TH-1
            do k=0,NZ_TH-1
              do i=0,NX_TH-1
                TH(i,k,j,ith-3)=tmp_th(i+1,j+1,k+1)
              end do
            end do
          end do

! Close dataset
          call h5sclose_f(filspace_id_th, error)
          call h5dclose_f(dset_id, error)
        end if

      end do

! Close the dataspace for the memory
      call h5sclose_f(memspace_id, error)
      call h5sclose_f(memspace_id_th, error)

      ! Close file
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      end subroutine ReadHDF5
