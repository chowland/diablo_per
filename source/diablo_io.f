C This file contains subroutines for inputting and outputting data in
C Diablo as well as all subroutines called directly from diablo.f
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INITIALIZE
      INCLUDE 'header'

      REAL    VERSION, CURRENT_VERSION
      logical RESET_TIME
      INTEGER I, J, K, N


      OPEN (11,file='input.dat',form='formatted',status='old')

C Read input file.
C   (Note - if you change the following section of code, update the
C    CURRENT_VERSION number to make obsolete previous input files!)

      CURRENT_VERSION=2.0
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) FLAVOR,   VERSION
      IF (VERSION .NE. CURRENT_VERSION) STOP 'Wrong input data format.'
      READ(11,*)
      READ(11,*) NU, LX, LY, LZ
      READ(11,*)
      READ(11,*) CREATE_NEW_FLOW
      READ(11,*)
      READ(11,*) N_TIME_STEPS, TIME_LIMIT, DELTA_T, RESET_TIME,
     &             VARIABLE_DT, CFL, UPDATE_DT
      READ(11,*)
      READ(11,*) VERBOSITY, SAVE_FLOW_INT, SAVE_STATS_INT, MOVIE
      READ(11,*)
      READ(11,*) NX_MOV, NX_MOV_TH, NY_MOV, NY_MOV_TH, NZ_MOV,
     &               NZ_MOV_TH
      READ(11,*)
! Read in the parameters for the N_TH scalars
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) CREATE_NEW_TH(N)
        READ(11,*)
        READ(11,*) RI_TAU(N), PR(N), REACTION(N)
      END DO

      LX=4.D0*ATAN(1.D0)*LX
      LY=4.D0*ATAN(1.D0)*LY
      LZ=4.D0*ATAN(1.D0)*LZ

C If we are using MPI, then Initialize the MPI Variables
      MPI_IO_NUM=''
      MPI_NUM=0
      CALL INIT_MPI
      IF (N_TH.ge.1) THEN
        CALL INIT_MPI_TH
      END IF

      IF (MOVIE) THEN
        CALL INIT_MOVIE
        if (rank.eq.0) write(*,*) 'Movie files created & initialized.'
      END IF
      CALL INIT_STATS
      if (rank.eq.0) write(*,*) 'Stats file created & initialized.'
      call init_mean
      if (rank.eq.0) write(*,*) 'H-mean file created & initialized.'
      call INIT_SPECTRA
      if (rank.eq.0) then
        call system('mkdir restart_files')
        write(*,*) 'Restart file directory created.'
      end if

C Initialize case-specific packages.
      CALL INPUT_PER
      CALL CREATE_GRID_PER
      CALL INIT_PER

C Initialize grid
      IF ((FLAVOR.NE.'Ensemble') .AND. (RANK.eq.0)) THEN
	      WRITE(6,*) 'Note that this code is distributed under the ',
     &           'GNU General Public License.'
        WRITE(6,*) 'No warranty is expressed or implied.'
        WRITE(6,*)
        write(*,*) 'Flavor: ',FLAVOR
        WRITE(6,*) 'Grid size: NX =',NX,', NY =',NY,', NZ =',NZ,'.'
        DO N=1,N_TH
          WRITE(6,*) 'Scalar number: ',N
          WRITE(6,*) '  Richardson number: ',RI_TAU(N)
          WRITE(6,*) '  Prandtl number: ',PR(N)
        END DO
	    END IF
      NXM=NX-1
      NYM=NY-1
      NZM=NZ-1
      NXM_TH=NX_TH-1
      NYM_TH=NY_TH-1
      NZM_TH=NZ_TH-1

C Initialize storage arrays.
      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NX+1
            U1(I,K,J)=0.
            U3(I,K,J)=0.
            U2(I,K,J)=0.
            P (I,K,J)=0.
            R1(I,K,J)=0.
            R2(I,K,J)=0.
            R3(I,K,J)=0.
            F1(I,K,J)=0.
            F2(I,K,J)=0.
            F3(I,K,J)=0.
          END DO
        END DO
      END DO
      DO J=0,NY+1
        DO K=0,NZ_S
          DO I=0,NX_S/2
            CU1(I,K,J)=0.
            CU3(I,K,J)=0.
            CU2(I,K,J)=0.
            CP (I,K,J)=0.
            CR1(I,K,J)=0.
            CR2(I,K,J)=0.
            CR3(I,K,J)=0.
            CF1(I,K,J)=0.
            CF2(I,K,J)=0.
            CF3(I,K,J)=0.
          END DO
        END DO
      END DO


C Initialize FFT package (includes defining the wavenumber vectors).
      write(*,*) 'Initializing FFT',RANK
      CALL INIT_FFT
      CALL INIT_FFT_TH
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      write(*,*) 'RANK, done init fft',RANK

C Initialize RKW3 parameters.
      H_BAR(1)=DELTA_T*(8.0/15.0)
      H_BAR(2)=DELTA_T*(2.0/15.0)
      H_BAR(3)=DELTA_T*(5.0/15.0)
      BETA_BAR(1)=1.0
      BETA_BAR(2)=25.0/8.0
      BETA_BAR(3)=9.0/4.0
      ZETA_BAR(1)=0.0
      ZETA_BAR(2)=-17.0/8.0
      ZETA_BAR(3)=-5.0/4.0

C Initialize values for reading of scalars
      NUM_READ_TH=0
      DO N=1,N_TH
        IF (CREATE_NEW_TH(N)) THEN
          NUM_READ_TH=NUM_READ_TH
        ELSE
          NUM_READ_TH=NUM_READ_TH + 1
          READ_TH_INDEX(NUM_READ_TH)=N
        END IF
      END DO
      CALL CREATE_TH_PER

C Create flow.
      IF (CREATE_NEW_FLOW) THEN
        CALL CREATE_FLOW_PER
        IF (FLAVOR.NE.'Ensemble') THEN
            write(*,*) 'A new flowfield has been created'
        END IF
        IF (FLAVOR.EQ.'Basic') THEN
            write(*,*) 'In DIABLO_IO...'
            CALL SAVE_STATS(.FALSE.)
            CALL SAVE_FLOW(.FALSE.)
        END IF

      ELSE
        CALL READ_FLOW

C Initialize flow.
        IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN
          PREVIOUS_TIME_STEP=0
          TIME_STEP=0
          TIME=0
        END IF

        CALL SAVE_STATS(.FALSE.)
        CALL POISSON_P_PER

      END IF

      IF (FLAVOR.eq.'CHEMOTAXIS' ) THEN
        EK=0.d0
        DO J=0,TNKY
          DO K=0,TNKZ_S
            DO I=0,NKX_S
              IF ((SQRT(KX2_S(I)+KY2(J)+KZ2_S(K)).le.2.5d0)
     &        .AND.((I+MOD(RANK,NP_S)*(NKX_S+1)).le.NKX)
     &        .AND.((K+INT(RANK/NP_S)*(TNKZ_S+1)).le.TNKZ)) THEN
                EK=EK+CU1(I,K,J)*CONJG(CU1(I,K,J))
     &               +CU2(I,K,J)*CONJG(CU2(I,K,J))
     &               +CU3(I,K,J)*CONJG(CU3(I,K,J))
              END IF
            END DO
          END DO
        END DO
        CALL MPI_ALLREDUCE(EK,EK_sum,1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
        EK0=EK_sum
        IF (RANK.EQ.0) WRITE(*,*) 'Initial EK0: ',EK0
      END IF

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      LOGICAL FINAL

      CALL SAVE_STATS_PER(FINAL)

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE READ_FLOW
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*60 FNAME
      INTEGER I, J, K, N

      complex*16 ctemp(200)

      FNAME='start.h5'
      IF (RANK.EQ.0) WRITE(6,*) 'Reading flow from ',FNAME
      
      call mpi_barrier(MPI_COMM_WORLD,ierror)
      call ReadHDF5(FNAME)

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_FLOW(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*60 FNAME
      INTEGER      I, J, K, N
      LOGICAL      FINAL, SAVE_PRESSURE

      SAVE_PRESSURE=.FALSE.
      if (FINAL) then
        FNAME='end.h5'
        SAVE_PRESSURE=.TRUE.
      else
!        if (N_TIME_STEPS/SAVE_FLOW_INT < 10) then
!          FNAME='out'//char(int(TIME_STEP/SAVE_FLOW_INT)+48)//'.h5'
        if (N_TIME_STEPS/SAVE_FLOW_INT < 100) then
          FNAME='out'//char(floor(TIME_STEP/SAVE_FLOW_INT/10.)+48)
     &           //char(mod(TIME_STEP/SAVE_FLOW_INT,10)+48)//'.h5'
        else if (N_TIME_STEPS/SAVE_FLOW_INT < 1000) then
          FNAME='out'//char(floor(TIME_STEP/SAVE_FLOW_INT/100.)+48)
     &           //char(floor(mod(TIME_STEP/SAVE_FLOW_INT,100)/10.)+48)
     &           //char(mod(TIME_STEP/SAVE_FLOW_INT,10)+48)//'.h5'
        end if
      end if

      call mpi_barrier(MPI_COMM_WORLD,ierror)
      call WriteHDF5('restart_files/'//FNAME,SAVE_PRESSURE)

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE END_RUN(FLAG)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      LOGICAL FLAG,FILE_EXISTS

      FLAG=.FALSE.
      ! Check for the time
      call WALL_TIME(END_TIME)
      if (END_TIME-START_TIME.gt.TIME_LIMIT) THEN
        IF (RANK.EQ.0) THEN
          write(*,*) ' STOP because of wall-time hit!'
        END IF
        FLAG=.TRUE.
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine wall_time(wt)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c
c     Return wall-clock time as seconds after Jan. 1, 2018.
c     Support for leap year is not included anymore.
c
c     By using a 'save' statement, the wall-time after the first
c     call to the subroutine could be computed, but that is not
c     intended with the present subroutine (e.g. the history file)
c
      implicit none

      real*8 wt
      integer val(8),i,shift,day

      integer mon(12,2)
      data mon /
     &     31,28,31,30,31,30,31,31,30,31,30,31,
     &     31,29,31,30,31,30,31,31,30,31,30,31/
c
c     Get current date and time
c     val(1) : year
c     val(2) : month
c     val(3) : day
c     val(4) : difference to GMT
c     val(5) : hour
c     val(6) : minute
c     val(7) : second
c     val(8) : 1/1000 second
c
      call date_and_time(values=val)
c
c     Determine leap year
c
      if (mod(val(1),4).eq.0) then
         if (mod(val(1),100).eq.0) then
            if (mod(val(1),400).eq.0) then
               shift=2
            else
               shift=1
            end if
         else
            shift=2
         end if
      else
         shift = 1
      end if
c
c     Construct day of the year
c
      day = val(3)-1
      do i=1,val(2)-1
         day=day+mon(i,shift)
      end do
c
c     And compute wall-clock time
c
      wt = (val(1)-2018)*365*86400+
     &     day*86400+val(5)*3600+val(6)*60+val(7)+dble(val(8)/1000.d0)

      end