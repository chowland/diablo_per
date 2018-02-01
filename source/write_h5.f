C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE WriteHDF5(FNAME,SAVE_PRESSURE)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      use hdf5

      INCLUDE 'header'

      CHARACTER*55 FNAME
      LOGICAL FINAL,SAVE_PRESSURE

      REAL*8 tmp(NX,NY_S+1,NZ_S+1)
      REAL*8 tmp_th(NX_TH,NY_S_TH+1,NZ_S_TH+1)
c
c     HDF5 ------------------------------------------------------
c
c     Dataset names
      character(len=10) :: dname

c     Identifiers
      integer(hid_t) :: file_id, dset_id
      integer(hid_t) :: filspace_id, memspace_id
      integer(hid_t) :: filspace_id_th, memspace_id_th

!     Identifiers
      integer(hid_t) :: selspace_id
      integer(hid_t) :: plist_id_w,plist_id_d,plist_id_d_th

c     Dimensions in the memory and in the file
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

      dimsm(1:3) = (/NX,NY_S+1,NZ_S+1/)
      dimsf(1:3) = (/NX,NY,NZ/)
      block = dimsm
      dimsm_th(1:3) = (/NX_TH,NY_S_TH+1,NZ_S_TH+1/)
      dimsf_th(1:3) = (/NX_TH,NY_TH,NZ_TH/)
      block_th = dimsm_th

      IF (RANK.EQ.0)
     &     WRITE(6,*) 'Writing flow to ',FNAME
      !write(6,*) 'RANKY,RANKZ: ',RANKY,RANKZ

      chunk_dims = (/NX,1,NZ_S+1/)
      chunk_dims_th = (/NX_TH,1,NZ_S_TH+1/)

! Stride and count for number of rows and columns in each dimension
      stride = 1
      count  = 1

! Offset determined by the rank of a processor
      offset_f = (/ 0, RANKY*(NY_S+1), RANKZ*(NZ_S+1) /)
      offset_th = (/ 0, RANKY*(NY_S_TH+1), RANKZ*(NZ_S_TH+1) /)
      offset_m(1:3)=0

! Initialize interface
      call h5open_f(error)

! Setup file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id_d, error)
      call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_world,
     +     mpi_info_null, error)

! Create the file collectively
      call h5fcreate_f(trim(FNAME), H5F_ACC_TRUNC_F,
     &                 file_id, error, access_prp = plist_id_d)
      call h5pclose_f(plist_id_d, error)

      adims=(/3,2/)
      arank=2
      call h5screate_simple_f(arank,adims,tspace, error)

      ! -----------------------------
      ! Resolution
      call h5acreate_f(file_id,'Resolution',H5T_STD_I32LE,tspace,
     &                 aid, error)
      tint(1)=NX
      tint(2)=NY
      tint(3)=NZ
      tint(4)=NX_TH
      tint(5)=NY_TH
      tint(6)=NZ_TH
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,tint,adims,error)
      call h5aclose_f(aid, error)
      ! -----------------------------
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

      ! -----------------------------
      ! Extra info
      adims=(/1,80/)
      tdim=80
      call h5tcopy_f(H5T_FORTRAN_S1,ttype,error)
      call h5tset_size_f(ttype,tdim,error)
      call h5screate_simple_f(arank,adims,tspace, error)

      call h5acreate_f(file_id,'Info',ttype,tspace,aid,
     &                 error)
      namnbuf=' (Put here what you want) '
      call h5awrite_f(aid,ttype,trim(namnbuf),adims,
     &                  error)
      call h5aclose_f(aid, error)
      call h5sclose_f(tspace,error)
      call h5tclose_f(ttype,error)
      ! -----------------------------

      call h5screate_f(H5S_SCALAR_F,tspace,error)

c$$$      ! -----------------------------
c$$$      ! HDF5-saving version
c$$$      call h5acreate_f(file_id,'Version',H5T_STD_I32LE,tspace,aid,
c$$$     &                 error)
c$$$      call h5awrite_f(aid,H5T_NATIVE_INTEGER,hver,adims,error)
c$$$      call h5aclose_f(aid, error)
c$$$      ! -----------------------------


      ! -----------------------------
      ! Time & Timestep
      call h5acreate_f(file_id,'Time',H5T_IEEE_F64LE,tspace,
     &                 aid, error)
      call h5awrite_f(aid,H5T_NATIVE_DOUBLE,TIME,adims,error)
      call h5aclose_f(aid, error)

      call h5acreate_f(file_id,'Timestep',H5T_STD_I32LE,tspace,
     &                 aid, error)
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,TIME_STEP,adims,error)
      call h5aclose_f(aid, error)
      ! ----------------------------

      call h5sclose_f(tspace,error)

      ! Convert to physical space
      call fft_xzy_mpi_to_physical(CU1,U1)
      call fft_xzy_mpi_to_physical(CU2,U2)
      call fft_xzy_mpi_to_physical(CU3,U3)
      do ith=1,N_TH
         CSTH1(:,:,:)=CTH(:,:,:,ith)
         call fft_xzy_mpi_to_physical_th(CSTH1,STH1)
         TH(:,:,:,ith)=STH1(:,:,:)
      end do

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
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id_w, error)
      call h5pset_dxpl_mpio_f(plist_id_w, H5FD_MPIO_COLLECTIVE_F,
     +        error)

      do ith=1,3+N_TH

        select case(ith)
        case (1)
          do j=0,NY_S
            do k=0,NZ_S
              do i=0,NXM
                tmp(i+1,j+1,k+1)=U1(i,k,j)
              end do
            end do
          end do
          dname="U"
        case (2)
          do j=0,NY_S
            do k=0,NZ_S
              do i=0,NXM
                tmp(i+1,j+1,k+1)=U2(i,k,j)
              end do
            end do
          end do
          dname="V"
        case (3)
          do j=0,NY_S
            do k=0,NZ_S
              do i=0,NXM
                tmp(i+1,j+1,k+1)=U3(i,k,j)
              end do
            end do
          end do
          dname="W"
        case (4:)
          do j=0,NY_S_TH
            do k=0,NZ_S_TH
              do i=0,NXM_TH
                tmp_th(i+1,j+1,k+1)=TH(i,k,j,ith-3)
              end do
            end do
          end do
          dname="TH"//CHAR(ith+45)
        end select
        if (ith.LE.3) then
! If saving velocity, use velocity grid
          call h5dcreate_f(file_id, trim(dname), H5T_IEEE_F64LE,
     +        filspace_id, dset_id, error, dcpl_id = plist_id_d)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
          call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
     &        offset_f, count, error, stride, block)

          call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block)

! Write the dataset collectively
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,
     +        tmp,
     +        dimsm, error, file_space_id = filspace_id,
     +        mem_space_id = memspace_id, xfer_prp = plist_id_w)

! Close dateset
          call h5dclose_f(dset_id, error)
        else
! If saving scalar, use TH grid
          call h5dcreate_f(file_id, trim(dname), H5T_IEEE_F64LE,
     +        filspace_id_th, dset_id, error, dcpl_id = plist_id_d_th)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
          call h5sselect_hyperslab_f (filspace_id_th, H5S_SELECT_SET_F,
     &        offset_th, count, error, stride, block_th)

          call h5sselect_hyperslab_f (memspace_id_th, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block_th)

! Write the dataset collectively
          call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,
     +        tmp_th,
     +        dimsm_th, error, file_space_id = filspace_id_th,
     +        mem_space_id = memspace_id_th, xfer_prp = plist_id_w)

! Close dateset
          call h5dclose_f(dset_id, error)
        end if
      end do

!     In the case of saving for the pressure as well
      if (SAVE_PRESSURE) then
         call fft_xzy_mpi_to_physical(CP,P)

         do j=0,NY_S
           do k=0,NZ_S
             do i=0,NXM
               tmp(i+1,j+1,k+1)=P(i,k,j)
             end do
           end do
         end do
         dname="P"

         call h5dcreate_f(file_id, trim(dname), H5T_IEEE_F64LE,
     +        filspace_id, dset_id, error, dcpl_id = plist_id_d)

!     Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
         call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
     &        offset_f, count, error, stride, block)

         call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block)

!     Write the dataset collectively
         call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,
     +        tmp,
     +        dimsm, error, file_space_id = filspace_id,
     +        mem_space_id = memspace_id, xfer_prp = plist_id_w)

!     Close dateset
         call h5dclose_f(dset_id, error)
         call fft_xzy_mpi_to_fourier(P,CP)
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

      call fft_xzy_mpi_to_fourier(U1,CU1)
      call fft_xzy_mpi_to_fourier(U2,CU2)
      call fft_xzy_mpi_to_fourier(U3,CU3)
      do ith=1,N_TH
         STH1(:,:,:)=TH(:,:,:,ith)
         call fft_xzy_mpi_to_fourier_th(STH1,CSTH1)
         CTH(:,:,:,ith)=CSTH1(:,:,:)
      end do

      ! call mpi_finalize(ierror)
      ! stop

      end subroutine WriteHDF5

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine time_string(cdt)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

c     Construct string in the format '19-DEC-2005 22:47:06'

      implicit none

      integer i

      integer val(8)
      character*20 cdt
      character*3 monc
      character*2 day, hour, minute, sec
      character*4 year

      call date_and_time(values=val)

      if (val(2).eq.1) then
         monc  = 'JAN'
      else if (val(2).eq.2) then
         monc  = 'FEB'
      else if (val(2).eq.3) then
         monc  = 'MAR'
      else if (val(2).eq.4) then
         monc  = 'APR'
      else if (val(2).eq.5) then
         monc  = 'MAY'
      else if (val(2).eq.6) then
         monc  = 'JUN'
      else if (val(2).eq.7) then
         monc  = 'JUL'
      else if (val(2).eq.8) then
         monc  = 'AUG'
      else if (val(2).eq.9) then
         monc  = 'SEP'
      else if (val(2).eq.10) then
         monc  = 'OCT'
      else if (val(2).eq.11) then
         monc  = 'NOV'
      else if (val(2).eq.12) then
         monc  = 'DEC'
      else
         monc  = 'XXX'
      end if

      write(day,'(i0.2)') val(3)
      write(year,'(i4)') val(1)
      write(hour,'(i0.2)') val(5)
      write(minute,'(i0.2)') val(6)
      write(sec,'(i0.2)') val(7)

      cdt = day // '-' // monc // '-' // year // ' ' //
     &      hour // ':' // minute // ':' // sec

      end subroutine time_string

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ReadHDF5(FNAME)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      use hdf5

      INCLUDE 'header'

      CHARACTER*55 FNAME
      LOGICAL FINAL, READ_PRESSURE

      REAL*8 tmp(NX,NY_S+1,NZ_S+1)
      REAL*8 tmp_th(NX_TH,NY_S_TH+1,NZ_S_TH+1)

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

      dimsm(1:3) = (/NX,NY_S+1,NZ_S+1/)
      block = dimsm
      dimsm_th(1:3) = (/NX_TH,NY_S_TH+1,NZ_S_TH+1/)
      block_th = dimsm_th

!     Stride and count for number of rows and columns in each dimension
      stride = 1
      count  = 1

!     Offset determined by the rank of a processor
      offset_f = (/ 0, RANKY*(NY_S+1), RANKZ*(NZ_S+1) /)
      offset_th = (/ 0, RANKY*(NY_S_TH+1), RANKZ*(NZ_S_TH+1) /)
      offset_m(1:3)=0

!     Initialize interface
      call h5open_f(error)

!     Setup file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id_d, error)
      call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_world,
     +     mpi_info_null, error)

!     Open the file
      call h5fopen_f(trim(FNAME), H5F_ACC_RDONLY_F,
     &                 file_id, error, access_prp = plist_id_d)
      call h5pclose_f(plist_id_d, error)

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
        if (RANK.eq.0) then
          write(*,*) ' Error. File and program have ',
     &        'different resolutions. '
          write(*,*) ' Program: ', NX,NY,NZ,NX_TH,NY_TH,NZ_TH
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
      call h5aread_f(aid,H5T_NATIVE_DOUBLE,TIME,adims,error)
      call h5aclose_f(aid, error)

      call h5aopen_by_name_f(file_id,'.','Timestep',aid,error)
      call h5aread_f(aid,H5T_NATIVE_INTEGER,TIME_STEP,adims,error)
      call h5aclose_f(aid, error)
!     -----------------------------

!     Create property list for collective dataset read
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id_w, error)
      call h5pset_dxpl_mpio_f(plist_id_w, H5FD_MPIO_COLLECTIVE_F,
     +        error)
!     Dataspaces in memory
      call h5screate_simple_f(rHDF5, dimsm, memspace_id, error)
      call h5screate_simple_f(rHDF5, dimsm_th, memspace_id_th, error)

      do ith=1,3+N_TH
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
!          write(*,*) 'RANK, RANKY, RANKZ, offset_f, block, count, '//
!     &      'stride: ',
!     &            RANK, RANKY, RANKZ, offset_f,block,count,stride
! Read the dataset collectively
          call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmp,
     +        dimsm, error, file_space_id = filspace_id,
     +        mem_space_id = memspace_id, xfer_prp = plist_id_w)
          select case(ith)
          case (1)
          do j=0,NY_S
            do k=0,NZ_S
              do i=0,NXM
                U1(i,k,j)=tmp(i+1,j+1,k+1)
               end do
            end do
          end do
          case (2)
          do j=0,NY_S
            do k=0,NZ_S
              do i=0,NXM
                U2(i,k,j)=tmp(i+1,j+1,k+1)
               end do
            end do
          end do
          case (3)
          do j=0,NY_S
            do k=0,NZ_S
              do i=0,NXM
                U3(i,k,j)=tmp(i+1,j+1,k+1)
               end do
            end do
          end do
          end select

! Close dataset
          call h5sclose_f(filspace_id, error)
          call h5dclose_f(dset_id, error)

        else if (.NOT.CREATE_NEW_TH(max(1,ith-3))) then
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
     +        mem_space_id = memspace_id_th, xfer_prp = plist_id_w)
          do j=0,NY_S_TH
            do k=0,NZ_S_TH
              do i=0,NXM_TH
                TH(i,k,j,ith-3)=tmp_th(i+1,j+1,k+1)
               end do
            end do
          end do

! Close dataset
          call h5sclose_f(filspace_id_th, error)
          call h5dclose_f(dset_id, error)
        end if

      end do

! Decide whether to compute the pressure or to read
      call h5lexists_f(file_id, 'P', READ_PRESSURE, error)
      if (READ_PRESSURE) then
        dname="P"
        call h5dopen_f(file_id,trim(dname),dset_id,error)
        call h5dget_space_f(dset_id,filspace_id,error)

! Select hyperslab in the file.
        call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
     &        offset_f, count, error, stride, block)
        call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block)

! Write the dataset collectively
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmp,
     +        dimsm, error, file_space_id = filspace_id,
     +        mem_space_id = memspace_id, xfer_prp = plist_id_w)

        do j=0,NY_S
          do k=0,NZ_S
            do i=0,NXM
              P(i,k,j)=tmp(i+1,j+1,k+1)
             end do
          end do
        end do
! Close dataset
        call h5sclose_f(filspace_id, error)
        call h5dclose_f(dset_id, error)
        call fft_xzy_mpi_to_fourier(P,CP)
      end if

! Close the dataspace for the memory
      call h5sclose_f(memspace_id, error)
      call h5sclose_f(memspace_id_th, error)

! Close the properties for the reading
      call h5pclose_f(plist_id_w, error)

      ! Close file
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      IF (VARIABLE_DT) THEN
        CALL COURANT_MPI
      END IF

      ! Convert to Fourier space
      call fft_xzy_mpi_to_fourier(U1,CU1)
      call fft_xzy_mpi_to_fourier(U2,CU2)
      call fft_xzy_mpi_to_fourier(U3,CU3)
      do ith=1,N_TH
        if (.NOT.CREATE_NEW_TH(ith)) then
         S1(:,:,:)=TH(:,:,:,ith)
         call fft_xzy_mpi_to_fourier_th(S1,CS1)
         CTH(:,:,:,ith)=CS1(:,:,:)
        end if
      end do

      end subroutine ReadHDF5

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_MOVIE
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      use hdf5

      INCLUDE 'header'

      CHARACTER*55 fname

c     HDF5 ------------------------------------------------------

c     Identifiers
      integer(hid_t) :: file_id, gid, plist_id_d, i

      integer :: arank, error, COMM_MPI, RANKY_MOV, RANKZ_MOV
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
          COMM_MPI=mpi_comm_y
        case(2)
          fname='movie_xz.h5'
          COMM_MPI=mpi_comm_z
        case(3)
          fname='movie_yz.h5'
          COMM_MPI=mpi_comm_world
        end select
        RANKY_MOV = NY_MOV/(NY_S+1)
        RANKZ_MOV = NZ_MOV/(NZ_S+1)
        if ( ((i.eq.1).and.(RANKZ.EQ.RANKZ_MOV)) .OR.
     &       ((i.eq.2).and.(RANKY.EQ.RANKY_MOV)) .OR.
     &        (i.eq.3) ) then
! Setup file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id_d, error)
      call h5pset_fapl_mpio_f(plist_id_d, COMM_MPI,
     +     mpi_info_null, error)

! Create the file collectively
      call h5fcreate_f(fname, H5F_ACC_TRUNC_F,
     &                 file_id, error, access_prp = plist_id_d)
      call h5pclose_f(plist_id_d, error)

! Write file attributes:

      ! ----------------------------
      ! Resolution
      adims=(/3,2/)
      arank=2
      call h5screate_simple_f(arank,adims,tspace, error)
      call h5acreate_f(file_id,'Resolution',H5T_STD_I32LE,tspace,
     &                 aid, error)
      tint(1)=NX
      tint(2)=NY
      tint(3)=NZ
      tint(4)=NX_TH
      tint(5)=NY_TH
      tint(6)=NZ_TH
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,tint,adims,error)
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
      call h5screate_f(H5S_SCALAR_F,tspace,error)
      call h5acreate_f(file_id,'Samples',H5T_STD_I32LE,tspace,aid,error)
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
      call h5aclose_f(aid,error)
      call h5sclose_f(tspace,error)

! Close file and interface
      call h5fclose_f(file_id,error)
      end if
      end do
      call h5close_f(error)

      end subroutine INIT_MOVIE

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE WriteHDF5_xyplane(IZ)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      use hdf5

      INCLUDE 'header'

      REAL*8 tmp(NX,NY_S+1)
      REAL*8 tmp_th(NX_TH,NY_S_TH+1)
c
c     HDF5 ------------------------------------------------------
c
c     Dataset names
      character(len=10) :: gname, dname

c     Identifiers
      integer(hid_t) :: file_id, dset_id
      integer(hid_t) :: filspace_id, memspace_id
      integer(hid_t) :: filspace_id_th, memspace_id_th

!     Identifiers
      integer(hid_t) :: gt_id, selspace_id
      integer(hid_t) :: plist_id_w,plist_id_d,plist_id_d_th

c     Dimensions in the memory and in the file
      integer(hsize_t), dimension(2) :: dimsm,dimsf
      integer(hsize_t), dimension(2) :: dimsm_th,dimsf_th

      integer(hsize_t), dimension(2) :: chunk_dims, count, offset_f
      integer(hsize_t), dimension(2) :: stride, block, offset_m
      integer(hsize_t), dimension(2) :: chunk_dims_th, offset_th
      integer(hsize_t), dimension(2) :: block_th

      integer :: rHDF5 = 2, arank, NSAMP, IZ, i, j
      logical flage

      integer(hsize_t),dimension(2) :: adims
      integer(size_t)               :: tdim
      integer(hid_t)                :: aid,tspace,ttype
      integer, dimension(1)         :: tint(6)
      character*80                  :: namnbuf
      character*20                  :: sttimec

      integer error, ith

! #### DEFINE WRITING PARAMETERS ####
      dimsm = (/NX,NY_S+1/)
      dimsf = (/NX,NY/)
      block = dimsm
      dimsm_th = (/NX_TH,NY_S_TH+1/)
      dimsf_th = (/NX_TH,NY_TH/)
      block_th = dimsm_th

      chunk_dims = (/NX,1/)
      chunk_dims_th = (/NX_TH,1/)

! Stride and count for number of rows and columns in each dimension
      stride = 1
      count  = 1

! Offset determined by the rank of a processor
      offset_f = (/ 0, RANKY*(NY_S+1) /)
      offset_th = (/ 0, RANKY*(NY_S_TH+1) /)
      offset_m = 0

! #### SET UP GROUP STRUCTURE ####

! Initialize interface
      call h5open_f(error)

! Setup file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id_d, error)
      call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_y,
     +     mpi_info_null, error)

! Open the file collectively and the relevant group
      call h5fopen_f('movie_xy.h5', H5F_ACC_RDWR_F,
     &                 file_id, error, access_prp = plist_id_d)
      call h5pclose_f(plist_id_d, error)

! Check and update the number of samples
      adims=1
      call h5aopen_f(file_id,'Samples',aid,error)
      call h5aread_f(aid,H5T_NATIVE_INTEGER,NSAMP,adims,error)
! Create group for this timestep
      write(gname,'(1i0.4)') NSAMP
      call h5gcreate_f(file_id,gname,gt_id,error)
      NSAMP=NSAMP+1
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,NSAMP,adims,error)
      call h5aclose_f(aid,error)

! Write Time, Timestep, missing coordinate attributes
      call h5screate_f(H5S_SCALAR_F,tspace,error)

      call h5acreate_f(gt_id,'Time',H5T_IEEE_F64LE,tspace,
     &                 aid, error)
      call h5awrite_f(aid,H5T_NATIVE_DOUBLE,TIME,adims,error)
      call h5aclose_f(aid, error)

      call h5acreate_f(gt_id,'Timestep',H5T_STD_I32LE,tspace,
     &                 aid, error)
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,TIME_STEP,adims,error)
      call h5aclose_f(aid, error)

      call h5acreate_f(gt_id,'z',H5T_IEEE_F64LE,tspace,aid,error)
      call h5awrite_f(aid,H5T_NATIVE_DOUBLE,GZ(IZ),adims,error)
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
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id_w, error)
      call h5pset_dxpl_mpio_f(plist_id_w, H5FD_MPIO_COLLECTIVE_F,
     +        error)
      if (RANKZ==NZ_MOV/(NZ_S+1)) then
      do ith=1,3+N_TH

      select case(ith)
      case (1)
        do i=0,NXM
          do j=0,NY_S
            tmp(i+1,j+1) = U1(i,IZ,j)
          end do
        end do
        dname="U"
      case (2)
        do i=0,NXM
          do j=0,NY_S
            tmp(i+1,j+1) = U2(i,IZ,j)
          end do
        end do
        dname="V"
      case (3)
        do i=0,NXM
          do j=0,NY_S
            tmp(i+1,j+1) = U3(i,IZ,j)
          end do
        end do
        dname="W"
      case (4:)
        do i=0,NXM_TH
          do j=0,NY_S_TH
            tmp_th(i+1,j+1) = TH(i,IZ,j,ith-3)
          end do
        end do
        dname="TH"//CHAR(ith+45)
      end select
!      if (ith.eq.1) write(*,*) 'RANK ',RANK,U1(:,IZ,NY_S)
      if (ith.LE.3) then
! If saving velocity, use velocity grid
        call h5dcreate_f(gt_id, trim(dname), H5T_IEEE_F64LE,
     +        filspace_id, dset_id, error, dcpl_id = plist_id_d)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
        call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
     &        offset_f, count, error, stride, block)

        call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block)

! Write the dataset collectively
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,
     +        tmp,
     +        dimsm, error, file_space_id = filspace_id,
     +        mem_space_id = memspace_id, xfer_prp = plist_id_w)

! Close dateset
        call h5dclose_f(dset_id, error)
      else
! If saving scalar, use TH grid
        call h5dcreate_f(gt_id, trim(dname), H5T_IEEE_F64LE,
     +        filspace_id_th, dset_id, error, dcpl_id = plist_id_d_th)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
        call h5sselect_hyperslab_f (filspace_id_th, H5S_SELECT_SET_F,
     &        offset_th, count, error, stride, block_th)

        call h5sselect_hyperslab_f (memspace_id_th, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block_th)

! Write the dataset collectively
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,
     +        tmp_th,
     +        dimsm_th, error, file_space_id = filspace_id_th,
     +        mem_space_id = memspace_id_th, xfer_prp = plist_id_w)

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

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE WriteHDF5_xzplane(IY)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      use hdf5

      INCLUDE 'header'


      REAL*8 tmp(NX,NZ_S+1)
      REAL*8 tmp_th(NX_TH,NZ_S_TH+1)
c
c     HDF5 ------------------------------------------------------
c
c     Dataset names
      character(len=10) :: gname, dname

c     Identifiers
      integer(hid_t) :: file_id, dset_id
      integer(hid_t) :: filspace_id, memspace_id
      integer(hid_t) :: filspace_id_th, memspace_id_th

!     Identifiers
      integer(hid_t) :: gt_id, selspace_id
      integer(hid_t) :: plist_id_w,plist_id_d,plist_id_d_th

c     Dimensions in the memory and in the file
      integer(hsize_t), dimension(2) :: dimsm,dimsf
      integer(hsize_t), dimension(2) :: dimsm_th,dimsf_th

      integer(hsize_t), dimension(2) :: chunk_dims, count, offset_f
      integer(hsize_t), dimension(2) :: stride, block, offset_m
      integer(hsize_t), dimension(2) :: chunk_dims_th, offset_th
      integer(hsize_t), dimension(2) :: block_th

      integer :: rHDF5 = 2, arank = 1, NSAMP, IY, i, k
      logical flage

      integer(hsize_t),dimension(1) :: adims
      integer(size_t)               :: tdim
      integer(hid_t)                :: aid,tspace,ttype
      integer, dimension(1)         :: tint(6)
      character*80                  :: namnbuf
      character*20                  :: sttimec

      integer error, ith

! #### DEFINE WRITING PARAMETERS ####
      dimsm = (/NX,NZ_S+1/)
      dimsf = (/NX,NZ/)
      block = dimsm
      dimsm_th = (/NX_TH,NZ_S_TH+1/)
      dimsf_th = (/NX_TH,NZ_TH/)
      block_th = dimsm_th

      chunk_dims = (/NX,NZ_S+1/)
      chunk_dims_th = (/NX_TH,NZ_S_TH+1/)

! Stride and count for number of rows and columns in each dimension
      stride = 1
      count  = 1

! Offset determined by the rank of a processor
      offset_f = (/ 0, RANKZ*(NZ_S+1) /)
      offset_th = (/ 0, RANKZ*(NZ_S_TH+1) /)
      offset_m = 0

! #### SET UP GROUP STRUCTURE ####

! Initialize interface
      call h5open_f(error)

! Setup file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id_d, error)
      call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_z,
     +     mpi_info_null, error)

! Open the file collectively and the relevant group
      call h5fopen_f('movie_xz.h5', H5F_ACC_RDWR_F,
     &                 file_id, error, access_prp = plist_id_d)
      call h5pclose_f(plist_id_d, error)

! Check and update the number of samples
      call h5aopen_f(file_id,'Samples',aid,error)
      call h5aread_f(aid,H5T_NATIVE_INTEGER,NSAMP,adims,error)
! Create group for this timestep
      write(gname,'(1i0.4)') NSAMP
      call h5gcreate_f(file_id,gname,gt_id,error)
      NSAMP=NSAMP+1
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,NSAMP,adims,error)
      call h5aclose_f(aid,error)

! Write Time, Timestep, missing coordinate attributes
      call h5screate_f(H5S_SCALAR_F,tspace,error)

      call h5acreate_f(gt_id,'Time',H5T_IEEE_F64LE,tspace,
     &                 aid, error)
      call h5awrite_f(aid,H5T_NATIVE_DOUBLE,TIME,adims,error)
      call h5aclose_f(aid, error)

      call h5acreate_f(gt_id,'Timestep',H5T_STD_I32LE,tspace,
     &                 aid, error)
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,TIME_STEP,adims,error)
      call h5aclose_f(aid, error)

      call h5acreate_f(gt_id,'y',H5T_IEEE_F64LE,tspace,aid,error)
      call h5awrite_f(aid,H5T_NATIVE_DOUBLE,GY(IY),adims,error)
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
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id_w, error)
      call h5pset_dxpl_mpio_f(plist_id_w, H5FD_MPIO_COLLECTIVE_F,
     +        error)
      if (RANKY==NY_MOV/(NY_S+1)) then
      do ith=1,3+N_TH

      select case(ith)
      case (1)
        do i=0,NXM
          do k=0,NZ_S
            tmp(i+1,k+1) = U1(i,k,IY)
          end do
        end do
        dname="U"
      case (2)
        do i=0,NXM
          do k=0,NZ_S
            tmp(i+1,k+1) = U2(i,k,IY)
          end do
        end do
        dname="V"
      case (3)
        do i=0,NXM
          do k=0,NZ_S
            tmp(i+1,k+1) = U3(i,k,IY)
          end do
        end do
        dname="W"
      case (4:)
        do i=0,NXM_TH
          do k=0,NZ_S_TH
            tmp_th(i+1,k+1) = TH(i,k,IY,ith-3)
          end do
        end do
        dname="TH"//CHAR(ith+45)
      end select
      if (ith.LE.3) then
! If saving velocity, use velocity grid
        call h5dcreate_f(gt_id, trim(dname), H5T_IEEE_F64LE,
     +        filspace_id, dset_id, error, dcpl_id = plist_id_d)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
        call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
     &        offset_f, count, error, stride, block)

        call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block)

! Write the dataset collectively
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,
     +        tmp,
     +        dimsm, error, file_space_id = filspace_id,
     +        mem_space_id = memspace_id, xfer_prp = plist_id_w)

! Close dateset
        call h5dclose_f(dset_id, error)
      else
! If saving scalar, use TH grid
        call h5dcreate_f(gt_id, trim(dname), H5T_IEEE_F64LE,
     +        filspace_id_th, dset_id, error, dcpl_id = plist_id_d_th)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
        call h5sselect_hyperslab_f (filspace_id_th, H5S_SELECT_SET_F,
     &        offset_th, count, error, stride, block_th)

        call h5sselect_hyperslab_f (memspace_id_th, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block_th)

! Write the dataset collectively
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,
     +        tmp_th,
     +        dimsm_th, error, file_space_id = filspace_id_th,
     +        mem_space_id = memspace_id_th, xfer_prp = plist_id_w)

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



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE WriteHDF5_yzplane(IX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      use hdf5

      INCLUDE 'header'

      REAL*8 tmp(NY_S+1,NZ_S+1)
      REAL*8 tmp_th(NY_S_TH+1,NZ_S_TH+1)
c
c     HDF5 ------------------------------------------------------
c
c     Dataset names
      character(len=10) :: gname, dname

c     Identifiers
      integer(hid_t) :: file_id, dset_id
      integer(hid_t) :: filspace_id, memspace_id
      integer(hid_t) :: filspace_id_th, memspace_id_th

!     Identifiers
      integer(hid_t) :: gt_id, selspace_id
      integer(hid_t) :: plist_id_w,plist_id_d,plist_id_d_th

c     Dimensions in the memory and in the file
      integer(hsize_t), dimension(2) :: dimsm,dimsf
      integer(hsize_t), dimension(2) :: dimsm_th,dimsf_th

      integer(hsize_t), dimension(2) :: chunk_dims, count, offset_f
      integer(hsize_t), dimension(2) :: stride, block, offset_m
      integer(hsize_t), dimension(2) :: chunk_dims_th, offset_th
      integer(hsize_t), dimension(2) :: block_th

      integer :: rHDF5 = 2, arank = 1, NSAMP, IX, j, k
      logical flage

      integer(hsize_t),dimension(1) :: adims
      integer(size_t)               :: tdim
      integer(hid_t)                :: aid,tspace,ttype
      integer, dimension(1)         :: tint(6)
      character*80                  :: namnbuf
      character*20                  :: sttimec

      integer error, ith

! #### DEFINE WRITING PARAMETERS ####
      dimsm = (/NY_S+1,NZ_S+1/)
      dimsf = (/NY,NZ/)
      block = dimsm
      dimsm_th = (/NY_S_TH+1,NZ_S_TH+1/)
      dimsf_th = (/NY_TH,NZ_TH/)
      block_th = dimsm_th

      chunk_dims = (/1,NZ_S+1/)
      chunk_dims_th = (/1,NZ_S_TH+1/)

! Stride and count for number of rows and columns in each dimension
      stride = 1
      count  = 1

! Offset determined by the rank of a processor
      offset_f = (/ RANKY*(NY_S+1), RANKZ*(NZ_S+1) /)
      offset_th = (/ RANKY*(NY_S_TH+1), RANKZ*(NZ_S_TH+1) /)
      offset_m = 0
! #### SET UP GROUP STRUCTURE ####

! Initialize interface
      call h5open_f(error)

! Setup file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id_d, error)
      call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_world,
     +     mpi_info_null, error)
! Open the file collectively and the relevant group
      call h5fopen_f('movie_yz.h5', H5F_ACC_RDWR_F,
     &                 file_id, error, access_prp = plist_id_d)
      call h5pclose_f(plist_id_d, error)

! Check and update the number of samples
      call h5aopen_f(file_id,'Samples',aid,error)
      call h5aread_f(aid,H5T_NATIVE_INTEGER,NSAMP,adims,error)
! Create group for this timestep
      write(gname,'(1i0.4)') NSAMP
      call h5gcreate_f(file_id,gname,gt_id,error)
      NSAMP=NSAMP+1
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,NSAMP,adims,error)
      call h5aclose_f(aid,error)
! Write Time, Timestep, missing coordinate attributes
      call h5screate_f(H5S_SCALAR_F,tspace,error)

      call h5acreate_f(gt_id,'Time',H5T_IEEE_F64LE,tspace,
     &                 aid, error)
      call h5awrite_f(aid,H5T_NATIVE_DOUBLE,TIME,adims,error)
      call h5aclose_f(aid, error)

      call h5acreate_f(gt_id,'Timestep',H5T_STD_I32LE,tspace,
     &                 aid, error)
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,TIME_STEP,adims,error)
      call h5aclose_f(aid, error)
      call h5acreate_f(gt_id,'x',H5T_IEEE_F64LE,tspace,aid,error)
      call h5awrite_f(aid,H5T_NATIVE_DOUBLE,GX(IX),adims,error)
      call h5aclose_f(aid,error)

      call h5sclose_f(tspace,error)
      ! ----------------------------
      call mpi_barrier(mpi_comm_world,ierr)
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
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id_w, error)
      call h5pset_dxpl_mpio_f(plist_id_w, H5FD_MPIO_COLLECTIVE_F,
     +        error)

      do ith=1,3+N_TH

      select case(ith)
      case (1)
        do j=0,NY_S
          do k=0,NZ_S
            tmp(j+1,k+1) = U1(IX,k,j)
          end do
        end do
        dname="U"
      case (2)
        do j=0,NY_S
          do k=0,NZ_S
            tmp(j+1,k+1) = U2(IX,k,j)
          end do
        end do
        dname="V"
      case (3)
        do j=0,NY_S
          do k=0,NZ_S
            tmp(j+1,k+1) = U3(IX,k,j)
          end do
        end do
        dname="W"
      case (4:)
        do j=0,NY_S_TH
          do k=0,NZ_S_TH
            tmp_th(j+1,k+1) = TH(IX,k,j,ith-3)
          end do
        end do
        dname="TH"//CHAR(ith+45)
      end select
      if (ith.LE.3) then
! If saving velocity, use velocity grid
        call h5dcreate_f(gt_id, trim(dname), H5T_IEEE_F64LE,
     +        filspace_id, dset_id, error, dcpl_id = plist_id_d)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
        call h5sselect_hyperslab_f (filspace_id, H5S_SELECT_SET_F,
     &        offset_f, count, error, stride, block)

        call h5sselect_hyperslab_f (memspace_id, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block)

! Write the dataset collectively
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,
     +        tmp,
     +        dimsm, error, file_space_id = filspace_id,
     +        mem_space_id = memspace_id, xfer_prp = plist_id_w)

! Close dateset
        call h5dclose_f(dset_id, error)
      else
! If saving scalar, use TH grid
        call h5dcreate_f(gt_id, trim(dname), H5T_IEEE_F64LE,
     +        filspace_id_th, dset_id, error, dcpl_id = plist_id_d_th)

! Select hyperslab in the file.
!     call h5dget_space_f(dsetur_id, selspace_id, error)
        call h5sselect_hyperslab_f (filspace_id_th, H5S_SELECT_SET_F,
     &        offset_th, count, error, stride, block_th)

        call h5sselect_hyperslab_f (memspace_id_th, H5S_SELECT_SET_F,
     +        offset_m, count, error, stride, block_th)

! Write the dataset collectively
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,
     +        tmp_th,
     +        dimsm_th, error, file_space_id = filspace_id_th,
     +        mem_space_id = memspace_id_th, xfer_prp = plist_id_w)

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


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_STATS
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      use hdf5

      INCLUDE 'header'

c     HDF5 ------------------------------------------------------

c     Identifiers
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
      call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_world,
     +     mpi_info_null, error)

! Create the file collectively
      call h5fcreate_f('stats.h5', H5F_ACC_TRUNC_F,
     &                 file_id, error, access_prp = plist_id_d)
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

      call h5acreate_f(file_id,'Date',ttype,tspace,aid,
     &                 error)
      call time_string(sttimec)
      call h5awrite_f(aid,ttype,trim(sttimec),adims,
     &                  error)
      call h5aclose_f(aid, error)
      call h5sclose_f(tspace,error)
      call h5tclose_f(ttype,error)
      ! -----------------------------

! Create the groups for each quantity (with attribute)
      call h5gcreate_f(file_id,'urms',gid,error)
      call h5screate_f(H5S_SCALAR_F,tspace,error)
      call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
      call h5aclose_f(aid,error)
      call h5sclose_f(tspace,error)
      call h5gclose_f(gid,error)

      call h5gcreate_f(file_id,'vrms',gid,error)
      call h5screate_f(H5S_SCALAR_F,tspace,error)
      call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
      call h5aclose_f(aid,error)
      call h5sclose_f(tspace,error)
      call h5gclose_f(gid,error)

      call h5gcreate_f(file_id,'wrms',gid,error)
      call h5screate_f(H5S_SCALAR_F,tspace,error)
      call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
      call h5aclose_f(aid,error)
      call h5sclose_f(tspace,error)
      call h5gclose_f(gid,error)

      call h5gcreate_f(file_id,'thth',gid,error)
      call h5screate_f(H5S_SCALAR_F,tspace,error)
      call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
      call h5aclose_f(aid,error)
      call h5sclose_f(tspace,error)
      call h5gclose_f(gid,error)

      call h5gcreate_f(file_id,'thvar',gid,error)
      call h5screate_f(H5S_SCALAR_F,tspace,error)
      call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
      call h5aclose_f(aid,error)
      call h5sclose_f(tspace,error)
      call h5gclose_f(gid,error)

      call h5gcreate_f(file_id,'thme',gid,error)
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

      call h5gcreate_f(file_id,'eta',gid,error)
      call h5screate_f(H5S_SCALAR_F,tspace,error)
      call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
      call h5aclose_f(aid,error)
      call h5sclose_f(tspace,error)
      call h5gclose_f(gid,error)

      call h5gcreate_f(file_id,'re_lambda',gid,error)
      call h5screate_f(H5S_SCALAR_F,tspace,error)
      call h5acreate_f(gid,'Samples',H5T_STD_I32LE,tspace,aid,error)
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,0,adims,error)
      call h5aclose_f(aid,error)
      call h5sclose_f(tspace,error)
      call h5gclose_f(gid,error)

! Close file and interface
      call h5fclose_f(file_id,error)
      call h5close_f(error)

      end subroutine INIT_STATS


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE WriteStatH5(gname,Stat)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      use hdf5

      INCLUDE 'header'

      REAL*8 Stat
c
c     HDF5 ------------------------------------------------------
c
c     Dataset names
      character(len=10) :: gname, dname

c     Identifiers
      integer(hid_t) :: file_id, dset_id
      integer(hid_t) :: filspace_id, memspace_id
      integer(hid_t) :: filspace_id_th, memspace_id_th

!     Identifiers
      integer(hid_t) :: gid, gt_id, selspace_id
!      integer(hid_t) :: plist_id_w,plist_id_d,plist_id_d_th

c     Dimensions in the memory and in the file
      integer(hsize_t), dimension(1) :: dims

      integer :: rHDF5 = 1, arank = 1, NSAMP
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

! Setup file access property list with parallel I/O access
!      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id_d, error)
!      call h5pset_fapl_mpio_f(plist_id_d, mpi_comm_world,
!     +     mpi_info_null, error)

! Open the file collectively and the relevant group
      call h5fopen_f('stats.h5', H5F_ACC_RDWR_F,file_id, error)
!      call h5pclose_f(plist_id_d, error)
      call h5gopen_f(file_id,'/'//gname,gid,error)

! Check and update the number of samples
      call h5aopen_f(gid,'Samples',aid,error)
      call h5aread_f(aid,H5T_NATIVE_INTEGER,NSAMP,adims,error)
! Create dataset name for this timestep
      write(dname,'(1i0.4)') NSAMP
      NSAMP=NSAMP+1
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,NSAMP,adims,error)
      call h5aclose_f(aid,error)

      call h5screate_f(H5S_SCALAR_F,tspace,error)

! Write statistic to file
      call h5dcreate_f(gid,dname,H5T_IEEE_F64LE,tspace,dset_id,error)
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,Stat,dims,error)

! Write Time, Timestep attributes
      call h5acreate_f(dset_id,'Time',H5T_IEEE_F64LE,tspace,
     &                 aid, error)
      call h5awrite_f(aid,H5T_NATIVE_DOUBLE,TIME,adims,error)
      call h5aclose_f(aid, error)

      call h5acreate_f(dset_id,'Timestep',H5T_STD_I32LE,tspace,
     &                 aid, error)
      call h5awrite_f(aid,H5T_NATIVE_INTEGER,TIME_STEP,adims,error)
      call h5aclose_f(aid, error)


      call h5dclose_f(dset_id,error)

      call h5sclose_f(tspace,error)
      ! ----------------------------

      ! Close groups
      call h5gclose_f(gid,error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      ! call mpi_finalize(ierror)
      ! stop

      end subroutine WriteStatH5
