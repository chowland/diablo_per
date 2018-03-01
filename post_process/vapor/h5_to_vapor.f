      program res_to_vapor
! This program converts a diablo restart file to a format readable by vapor
! A specified number of wavelet transforms are performed to produce
! compressed versions of the data
! Note, as of version 1.1.2, Vapor cannot handle stretched grids

! grid_def should contain the grid size and number of scalars
      include 'header'

      character*85 FNAME
      logical CREATE_NEW_FLOW

! Define Variables
      integer i,j,k,n
      real buffer3(0:NX-1,0:NZ-1,0:NY-1)
! The sizes of these strings may need to be changed
      character*15 size_str
      character*55 len_str,run_dir
      character*4 NUM_TIMES

! **** User Input *****
! Variables for movie
      integer FIRST_OUT, LAST_OUT, OUT_NUM, loop_index
      FIRST_OUT=0
      LAST_OUT=60
      NUM_TIMES='61' ! Remember to add 1 for start.h5 here if not CREATE_NEW_FLOW

! Number of periodic directions used in the simulation
      NUM_PER_DIR=3
      CREATE_NEW_FLOW=.TRUE.
! This string should contain the size of the buffer array
      size_str='128x128x128'
      len_str='0.0:0.0:0.0:6.28:6.28:6.28'
      LX=6.28
      LY=6.28
      LZ=6.28

      run_dir='/local/scratch/public/cjh225/GM9_init_forced' ! Run directory

!      RI_TAU(1)=1.0d0

      NXM=NX-1
      NYM=NY-1
      NZM=NZ-1

      call init_fft

         DO I=0,NXM
           GX(I)=(DBLE(I)*LX)/DBLE(NX)
         END DO
         DO J=0,NYM
           GY(J)=(DBLE(J)*LY)/DBLE(NY)
         END DO
         DO K=0,NZM
           GZ(K)=(DBLE(K)*LZ)/DBLE(NZ)
         END DO

      write(*,*) 'creating header...'

! Now we are ready to create a vdf (header) file
      if (N_TH.eq.0) then
        call SYSTEM('vdfcreate -dimension '
     &//size_str//' -extents '//len_str//' -periodic 1:1:1
     &    -numts '//NUM_TIMES//' -level 3 -vars3d U1:U2:U3
     &     vapor.vdf')
      else if (N_TH.eq.1) then
        call SYSTEM('vdfcreate -dimension '
     &//size_str//' -extents '//len_str//' -periodic 1:1:1
     &    -numts '//NUM_TIMES//' -level 3 -vars3d U1:U2:U3:TH1
     &      vapor.vdf')
      else if (N_TH.eq.2) then
        call SYSTEM('vdfcreate -dimension '
     &//size_str//' -extents '//len_str//' -periodic 1:1:1
     &    -numts '//NUM_TIMES//' -level 3 -vars3d U1:U2:U3:TH1:TH2
     &      vapor.vdf')
      else if (N_TH.eq.3) then
        call SYSTEM('vdfcreate -dimension '
     &//size_str//' -extents '//len_str//' -periodic 1:1:1
     &    -numts '//NUM_TIMES//' -level 3 -vars3d U1:U2:U3:TH1:TH2:TH3
     &      vapor.vdf')
      end if

      if (.NOT. CREATE_NEW_FLOW) then
        FNAME=trim(run_dir)//'/start.h5'
        write(*,*) 'reading flow from ',FNAME
        call ReadHDF5(FNAME)

        do i=0,NXM
          do j=0,NYM
            do k=0,NZM
              buffer3(i,k,j)=real(U1(i,k,j))
            end do
          end do
        end do
        open (22,file='temp.raw',form='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=NX*NY*NZ)
        write(22,REC=1) buffer3
        close(22)
        write(*,*) 'raw2vdf for U1...'
        call SYSTEM('raw2vdf -ts 0 -varname U1 vapor.vdf temp.raw')

        do i=0,NXM
          do j=0,NYM
            do k=0,NZM
              buffer3(i,k,j)=real(U2(i,k,j))
            end do
          end do
        end do
        open (22,file='temp.raw',form='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=NX*NY*NZ)
        write(22,REC=1) buffer3
        close(22)
        write(*,*) 'raw2vdf for U2...'
        call SYSTEM('raw2vdf -ts 0 -varname U2 vapor.vdf temp.raw')

        do i=0,NXM
          do j=0,NYM
            do k=0,NZM
              buffer3(i,k,j)=real(U3(i,k,j))
            end do
          end do
        end do
        open (22,file='temp.raw',form='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=NX*NY*NZ)
        write(22,REC=1) buffer3
        close(22)
        write(*,*) 'raw2vdf for U3...'
        call SYSTEM('raw2vdf -ts 0 -varname U3 vapor.vdf temp.raw')

        if (N_TH.ge.1) then
          do i=0,NXM
            do j=0,NYM
              do k=0,NZM
                buffer3(i,k,j)=real(TH(i,k,j,1))
              end do
            end do
          end do
          open (22,file='temp.raw',form='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=NX*NY*NZ)
          write(22,REC=1) buffer3
          close(22)
          write(*,*) 'raw2vdf for TH1...'
          call SYSTEM('raw2vdf -ts 0 -varname TH1 vapor.vdf temp.raw')
        end if

        if (N_TH.ge.2) then
          do i=0,NXM
            do j=0,NYM
              do k=0,NZM
                buffer3(i,k,j)=real(TH(i,k,j,2))
              end do
            end do
          end do
          open (22,file='temp.raw',form='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=NX*NY*NZ)
          write(22,REC=1) buffer3
          close(22)
          write(*,*) 'raw2vdf for TH2...'
          call SYSTEM('raw2vdf -ts 0 -varname TH2 vapor.vdf temp.raw')
        end if

        if (N_TH.ge.3) then
          do i=0,NXM
            do j=0,NYM
              do k=0,NZM
                buffer3(i,k,j)=real(TH(i,k,j,3))
              end do
            end do
          end do
          open (22,file='temp.raw',form='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=NX*NY*NZ)
          write(22,REC=1) buffer3
          close(22)
          write(*,*) 'raw2vdf for TH3...'
          call SYSTEM('raw2vdf -ts 0 -varname TH3 vapor.vdf temp.raw')
        end if

      end if

      DO OUT_NUM=FIRST_OUT,LAST_OUT
        if (CREATE_NEW_FLOW) then
          loop_index=OUT_NUM-FIRST_OUT
        else
          loop_index=OUT_NUM-FIRST_OUT+1
        end if
        if (LAST_OUT.lt.10) then
          FNAME=trim(run_dir)//'/restart_files/out'
     &        //char(mod(OUT_NUM,10)+48)//'.h5'
        else if (LAST_OUT.lt.100) then
          FNAME=trim(run_dir)//'/restart_files/out'
     &        //char(mod(OUT_NUM,100)/10+48)
     &        //char(mod(OUT_NUM,10)+48)//'.h5'
        else if (LAST_OUT.lt.1000) then
          FNAME=trim(run_dir)//'/restart_files/out'
     &        //char(mod(OUT_NUM,1000)/100+48)
     &        //char(mod(OUT_NUM,100)/10+48)
     &        //char(mod(OUT_NUM,10)+48)//'.h5'
        else
          stop 'Error: too many output files'
        end if
        write(*,*) 'reading flow from ',FNAME
        call ReadHDF5(FNAME)

        do i=0,NXM
          do j=0,NYM
            do k=0,NZM
              buffer3(i,k,j)=real(U1(i,k,j))
            end do
          end do
        end do
        open (22,file='temp.raw',form='UNFORMATTED',
     &    ACCESS='DIRECT',RECL=NX*NY*NZ)
        write(22,REC=1) buffer3
        close(22)
        write(*,*) 'raw2vdf for U1...'
        if (loop_index.lt.10) then
          call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,10)+48)
     &      //' -varname U1 vapor.vdf temp.raw')
        else if (loop_index.lt.100) then
          call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,100)/10+48)
     &      //char(mod(loop_index,10)+48)
     &      //' -varname U1 vapor.vdf temp.raw')
        else if (loop_index.lt.1000) then
          call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,1000)/100+48)
     &      //char(mod(loop_index,100)/10+48)
     &      //char(mod(loop_index,10)+48)
     &      //' -varname U1 vapor.vdf temp.raw')
        end if

        do i=0,NXM
          do j=0,NYM
            do k=0,NZM
              buffer3(i,k,j)=real(U2(i,k,j))
            end do
          end do
        end do
        open (22,file='temp.raw',form='UNFORMATTED',
     &    ACCESS='DIRECT',RECL=NX*NY*NZ)
        write(22,REC=1) buffer3
        close(22)
        write(*,*) 'raw2vdf for U2...'
        if (loop_index.lt.10) then
          call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,10)+48)
     &      //' -varname U2 vapor.vdf temp.raw')
        else if (loop_index.lt.100) then
          call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,100)/10+48)
     &      //char(mod(loop_index,10)+48)
     &      //' -varname U2 vapor.vdf temp.raw')
        else if (loop_index.lt.1000) then
          call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,1000)/100+48)
     &      //char(mod(loop_index,100)/10+48)
     &      //char(mod(loop_index,10)+48)
     &      //' -varname U2 vapor.vdf temp.raw')
        end if

        do i=0,NXM
          do j=0,NYM
            do k=0,NZM
              buffer3(i,k,j)=real(U3(i,k,j))
            end do
          end do
        end do
        open (22,file='temp.raw',form='UNFORMATTED',
     &    ACCESS='DIRECT',RECL=NX*NY*NZ)
        write(22,REC=1) buffer3
        close(22)
        write(*,*) 'raw2vdf for U3...'
        if (loop_index.lt.10) then
          call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,10)+48)
     &      //' -varname U3 vapor.vdf temp.raw')
        else if (loop_index.lt.100) then
          call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,100)/10+48)
     &      //char(mod(loop_index,10)+48)
     &      //' -varname U3 vapor.vdf temp.raw')
        else if (loop_index.lt.1000) then
          call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,1000)/100+48)
     &      //char(mod(loop_index,100)/10+48)
     &      //char(mod(loop_index,10)+48)
     &      //' -varname U3 vapor.vdf temp.raw')
        end if

        if (N_TH.ge.1) then
          do i=0,NXM
            do j=0,NYM
              do k=0,NZM
                buffer3(i,k,j)=real(TH(i,k,j,1))
              end do
            end do
          end do
          open (22,file='temp.raw',form='UNFORMATTED',
     &    ACCESS='DIRECT',RECL=NX*NY*NZ)
          write(22,REC=1) buffer3
          close(22)
          write(*,*) 'raw2vdf for TH1...'
          if (loop_index.lt.10) then
            call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,10)+48)
     &      //' -varname TH1 vapor.vdf temp.raw')
          else if (loop_index.lt.100) then
            call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,100)/10+48)
     &      //char(mod(loop_index,10)+48)
     &      //' -varname TH1 vapor.vdf temp.raw')
          else if (loop_index.lt.1000) then
            call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,1000)/100+48)
     &      //char(mod(loop_index,100)/10+48)
     &      //char(mod(loop_index,10)+48)
     &      //' -varname TH1 vapor.vdf temp.raw')
          end if
        end if

        if (N_TH.ge.2) then
          do i=0,NXM
            do j=0,NYM
              do k=0,NZM
                buffer3(i,k,j)=real(TH(i,k,j,2))
              end do
            end do
          end do
          open (22,file='temp.raw',form='UNFORMATTED',
     &    ACCESS='DIRECT',RECL=NX*NY*NZ)
          write(22,REC=1) buffer3
          close(22)
          write(*,*) 'raw2vdf for TH2...'
          if (loop_index.lt.10) then
            call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,10)+48)
     &      //' -varname TH2 vapor.vdf temp.raw')
          else if (loop_index.lt.100) then
            call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,100)/10+48)
     &      //char(mod(loop_index,10)+48)
     &      //' -varname TH2 vapor.vdf temp.raw')
          else if (loop_index.lt.1000) then
            call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,1000)/100+48)
     &      //char(mod(loop_index,100)/10+48)
     &      //char(mod(loop_index,10)+48)
     &      //' -varname TH2 vapor.vdf temp.raw')
          end if
        end if

        if (N_TH.ge.3) then
          do i=0,NXM
            do j=0,NYM
              do k=0,NZM
                buffer3(i,k,j)=real(TH(i,k,j,3))
              end do
            end do
          end do
          open (22,file='temp.raw',form='UNFORMATTED',
     &    ACCESS='DIRECT',RECL=NX*NY*NZ)
          write(22,REC=1) buffer3
          close(22)
          write(*,*) 'raw2vdf for TH3...'
          if (loop_index.lt.10) then
            call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,10)+48)
     &      //' -varname TH3 vapor.vdf temp.raw')
          else if (loop_index.lt.100) then
            call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,100)/10+48)
     &      //char(mod(loop_index,10)+48)
     &      //' -varname TH3 vapor.vdf temp.raw')
          else if (loop_index.lt.1000) then
            call SYSTEM('raw2vdf -ts '
     &      //char(mod(loop_index,1000)/100+48)
     &      //char(mod(loop_index,100)/10+48)
     &      //char(mod(loop_index,10)+48)
     &      //' -varname TH3 vapor.vdf temp.raw')
          end if
        end if

      end do

!      call SYSTEM('rm -f temp.raw')

!     We are done
      stop
      end
