C This file contains subroutines for inputting and outputting data in
C Diablo as well as all subroutines called directly from diablo.f
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INITIALIZE
      INCLUDE 'mpif.h'
      INCLUDE 'header'
      INCLUDE 'header_mpi'

      REAL    VERSION, CURRENT_VERSION
      logical RESET_TIME
      INTEGER I, J, K, N


      OPEN (11,file='input.dat',form='formatted',status='old')      

C Read input file.
C   (Note - if you change the following section of code, update the
C    CURRENT_VERSION number to make obsolete previous input files!)

      CURRENT_VERSION=1.0
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) FLAVOR,   VERSION
      IF (VERSION .NE. CURRENT_VERSION) STOP 'Wrong input data format.'
      READ(11,*)
      READ(11,*) USE_MPI,    LES
      READ(11,*)
      READ(11,*) NU, LX, LY, LZ
      READ(11,*)
      READ(11,*) NUM_PER_DIR, CREATE_NEW_FLOW
      READ(11,*)
      READ(11,*) N_TIME_STEPS, DELTA_T, RESET_TIME, VARIABLE_DT, CFL
     &            , UPDATE_DT
      READ(11,*)
      READ(11,*) VERBOSITY, SAVE_FLOW_INT, SAVE_STATS_INT, MOVIE
      READ(11,*)
! Read in the parameters for the N_TH scalars
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) CREATE_NEW_TH(N)
        READ(11,*)
        READ(11,*) FILTER_TH(N), FILTER_INT(N)
        READ(11,*)
        READ(11,*) RI_TAU(N), PR(N), REACTION(N)
      END DO

C If we are using MPI, then Initialize the MPI Variables
      MPI_IO_NUM=''
      MPI_NUM=0
      IF (USE_MPI) THEN
        CALL INIT_MPI
        IF (N_TH.ge.1) THEN
          CALL INIT_MPI_TH
        END IF
      END IF


C Initialize case-specific packages.
      IF (NUM_PER_DIR.EQ.3) THEN
        CALL INPUT_PER
        CALL CREATE_GRID_PER
        CALL INIT_PER
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        CALL INPUT_CHAN
        CALL CREATE_GRID_CHAN
        IF (USE_MPI) THEN
          CALL INIT_CHAN_MPI
        ELSE 
          CALL INIT_CHAN
        END IF
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        CALL INPUT_DUCT
        CALL CREATE_GRID_DUCT
        CALL INIT_DUCT
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        CALL INPUT_CAV 
        CALL CREATE_GRID_CAV
        CALL INIT_CAV
      END IF

C Initialize grid
      IF (FLAVOR.NE.'Ensemble') THEN
	WRITE(6,*) 'Note that this code is distributed under the ',
     *           'GNU General Public License.'
      WRITE(6,*) 'No warranty is expressed or implied.'
      WRITE(6,*)
      write(*,*) 'Flavor: ',FLAVOR
      WRITE(6,*) 'Grid size: NX =',NX,', NY =',NY,', NZ =',NZ,'.'
      DO N=1,N_TH
        WRITE(6,*) 'Scalar number: ',N
        WRITE(6,*) '  Richardson number: ',RI_TAU(N)
        WRITE(6,*) '  Prandlt number: ',PR(N)
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
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
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
      IF (NUM_PER_DIR.EQ.2) THEN
        CALL CREATE_TH_CHAN 
      ELSE IF (NUM_PER_DIR.EQ.3) THEN
        CALL CREATE_TH_PER
      END IF 
	
C Create flow.
      IF (CREATE_NEW_FLOW) THEN
        IF (RANK.eq.0) write(*,*) 'Creating flow...'
        IF (NUM_PER_DIR.EQ.3) THEN
          CALL CREATE_FLOW_PER
        ELSEIF (NUM_PER_DIR.EQ.2) THEN
          CALL CREATE_FLOW_CHAN
        ELSEIF (NUM_PER_DIR.EQ.1) THEN
          CALL CREATE_FLOW_DUCT
        ELSEIF (NUM_PER_DIR.EQ.0) THEN
          CALL CREATE_FLOW_CAV
        END IF
        IF (FLAVOR.NE.'Ensemble') THEN
            write(*,*) 'A new flowfield has been created'
        END IF
        IF (FLAVOR.EQ.'Basic') THEN
            write(*,*) 'In DIABLO_IO...'
            CALL SAVE_STATS(.FALSE.)
            CALL SAVE_FLOW(.FALSE.)
        END IF

      ELSE
        IF (RANK.eq.0) write(*,*) 'Reading flow...'
        CALL READ_FLOW

C Temporary...
!        CALL SAVE_FLOW(.FALSE.)
 
C Initialize flow.
      IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN
        PREVIOUS_TIME_STEP=0
        TIME_STEP=0
        TIME=0
      END IF

        CALL SAVE_STATS(.FALSE.)
        IF (NUM_PER_DIR.EQ.3) THEN
          CALL POISSON_P_PER
        ELSEIF (NUM_PER_DIR.EQ.2) THEN
          CALL POISSON_P_CHAN
        ELSEIF (NUM_PER_DIR.EQ.1) THEN
          CALL POISSON_P_DUCT
        ELSEIF (NUM_PER_DIR.EQ.0) THEN
          CALL POISSON_P_CAV
        END IF

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

      IF (NUM_PER_DIR.EQ.3) THEN
        CALL SAVE_STATS_PER(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        CALL SAVE_STATS_CHAN(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        CALL SAVE_STATS_DUCT(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        CALL SAVE_STATS_CAV(FINAL)          
      END IF

      IF (FINAL) THEN
        IF (NUM_PER_DIR.EQ.3) THEN
          CALL VIS_FLOW_PER         
        ELSEIF (NUM_PER_DIR.EQ.2) THEN
          CALL VIS_FLOW_CHAN         
        ELSEIF (NUM_PER_DIR.EQ.1) THEN
          CALL VIS_FLOW_DUCT          
        ELSEIF (NUM_PER_DIR.EQ.0) THEN
          CALL VIS_FLOW_CAV         
        END IF
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE READ_FLOW
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'mpif.h'
      INCLUDE 'header'
      include 'header_mpi'
      CHARACTER*60 FNAME
      CHARACTER*60 FNAME_P
      CHARACTER*60 FNAME_TH(N_TH)
      INTEGER I, J, K, N

      complex*16 ctemp(200)

           FNAME='diablo.start'
           FNAME_P='diablo_p.start'
           DO N=1,N_TH
             FNAME_TH(N)='diablo_th'
     &          //CHAR(MOD(N,100)/10+48)
     &          //CHAR(MOD(N,10)+48) // '.start'
           END DO

      WRITE(6,*)   'Reading flow from ',FNAME

      IF (.NOT.USE_MPI) THEN      
        OPEN(UNIT=10,FILE=FNAME,STATUS="OLD",FORM="UNFORMATTED")
        READ (10) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP
        READ (10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY)
! Check to make sure that the array dimensions match
      NKX_T=NX_T/3
      TNKZ_T=2*(NZ_T/3) 
      IF ((NX .NE. NX_T) .OR. (NY .NE. NY_T) .OR. (NZ .NE. NZ_T)) THEN
           write(*,*) 'NX,NY,NZ: ',NX,NY,NZ
           write(*,*) 'NX_T,NY_T,NZ_T: ',NX_T,NY_T,NZ_T
           STOP 'Error: old flowfield wrong dimensions. '
      END IF
      IF (NUM_PER_DIR .NE. NUM_PER_DIR_T)
     *     STOP 'Error: old flowfield wrong NUM_PER_DIR. '

      IF (NUM_PER_DIR.EQ.3) THEN
        DO N=1,NUM_READ_TH
! Specify in input.dat which scalars are to be read
          OPEN(UNIT=11,FILE=FNAME_TH(READ_TH_INDEX(N)),STATUS="OLD"
     &           ,FORM="UNFORMATTED")
          READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP
          READ (11) (((CTH(I,K,J,READ_TH_INDEX(N))
     &           ,I=0,NKX),K=0,TNKZ),J=0,TNKY)
         CLOSE(11)
        END DO
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        if ((NX.ne.NX_T).or.(NY.ne.NY_T).or.(NZ.ne.NZ_T)) then
! Interpolation....
        write(*,*) 'NX,NX_T,NKX ',NX,NX_t,NX_T/3
        write(*,*) 'NY,NY_T,TNKZ',NY,NY_T,2*(NZ_T/3)
        write(*,*) 'NZ,NZ_T',NZ,NZ_T
        write(*,*) 'interpolating'
        READ (10) (((CU1(I,K,J),I=0,NKX_T),K=0,TNKZ_T),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX_T),K=0,TNKZ_T),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX_T),K=0,TNKZ_T),J=1,NY)
        DO N=1,NUM_READ_TH
! Specify in input.dat which scalars are to be read
          OPEN(UNIT=11,FILE=FNAME_TH(READ_TH_INDEX(N)),STATUS="OLD"
     &           ,FORM="UNFORMATTED")
          READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP
          READ (11) (((CTH(I,K,J,READ_TH_INDEX(N))
     &           ,I=0,NKX_T),K=0,TNKZ_T),J=1,NY)
         CLOSE(11)
        END DO
        do j=1,NY
        do k=1,NZ_T/3
        do i=0,NKX
          CU1(I,TNKZ+1-k,J)=CU1(I,TNKZ_T+1-k,J)
          CU2(I,TNKZ+1-k,J)=CU2(I,TNKZ_T+1-k,J)
          CU3(I,TNKZ+1-k,J)=CU3(I,TNKZ_T+1-k,J)
          do N=1,NUM_READ_TH
            CTH(I,TNKZ+1-k,J,N)=CTH(I,TNKZ_T+1-k,J,N)
          end do
        end do
        end do
        end do
        do j=1,NY
        do k=NZ_T/3,TNKZ-NZ_T/3
        do i=0,NKX
          CU1(I,K,J)=0.d0
          CU2(I,K,J)=0.d0
          CU3(I,K,J)=0.d0
          do N=1,NUM_READ_TH
            CTH(I,K,J,N)=0.d0
          end do
         end do
         end do
         end do
        else
! No interpolation
        READ (10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY)
        DO N=1,NUM_READ_TH
! Specify in input.dat which scalars are to be read
          OPEN(UNIT=11,FILE=FNAME_TH(READ_TH_INDEX(N)),STATUS="OLD"
     &           ,FORM="UNFORMATTED")
          READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP
          READ (11) (((CTH(I,K,J,READ_TH_INDEX(N))
     &           ,I=0,NKX),K=0,TNKZ),J=1,NY)
         CLOSE(11)
        END DO
        end if
! Done if NUM_PER_DIR.EQ.2
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        READ (10) (((CU1(I,K,J),I=0,NKX),K=1,NZ),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX),K=1,NZ),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX),K=2,NZ),J=1,NY)
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        READ (10) (((U1(I,K,J),I=2,NX),K=1,NZ),J=1,NY),
     *            (((U2(I,K,J),I=1,NX),K=1,NZ),J=2,NY),
     *            (((U3(I,K,J),I=1,NX),K=2,NZ),J=1,NY)
      END IF
      CLOSE(10)
      CLOSE(11)
      ELSE

        CALL MPI_IO_READ(FNAME,FNAME_TH,FNAME_P)

!        CALL SAVE_STATS(.FALSE.)

      END IF


C Apply initial boundary conditions, set ghost cells
!      IF (USE_MPI) THEN
!        call APPLY_BC_VEL_MPI
!      ELSE
!        call APPLY_BC_VEL_LOWER
!        call APPLY_BC_VEL_UPPER
!      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_FLOW(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'mpif.h'
      INCLUDE 'header'
      include 'header_mpi'
      CHARACTER*60 FNAME
      CHARACTER*60 FNAME_P
      CHARACTER*60 FNAME_TH(N_TH)
      INTEGER      I, J, K, N
      LOGICAL      FINAL

      IF (FINAL) THEN
        FNAME='diablo.end'
        FNAME_P='diablo_p.end'
        DO N=1,N_TH
          FNAME_TH(N)='diablo_th'
     &       //CHAR(MOD(N,100)/10+48)
     &       //CHAR(MOD(N,10)+48) // '.end'
        END DO
      ELSE
         FNAME='./restart_files/diablo.saved.'
     &        //CHAR(MOD(TIME_STEP,100000)/10000+48)
     &        //CHAR(MOD(TIME_STEP,10000)/1000+48)
     &        //CHAR(MOD(TIME_STEP,1000)/100+48)
     &        //CHAR(MOD(TIME_STEP,100)/10+48)
     &        //CHAR(MOD(TIME_STEP,10)+48)
         FNAME_P='./restart_files/diablo_p.saved.'
     &        //CHAR(MOD(TIME_STEP,100000)/10000+48)
     &        //CHAR(MOD(TIME_STEP,10000)/1000+48)
     &        //CHAR(MOD(TIME_STEP,1000)/100+48)
     &        //CHAR(MOD(TIME_STEP,100)/10+48)
     &        //CHAR(MOD(TIME_STEP,10)+48)
           DO N=1,N_TH
           FNAME_TH(N)='./restart_files/diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) // '.saved.'
     &        //CHAR(MOD(TIME_STEP,100000)/10000+48)
     &        //CHAR(MOD(TIME_STEP,10000)/1000+48)
     &        //CHAR(MOD(TIME_STEP,1000)/100+48)
     &        //CHAR(MOD(TIME_STEP,100)/10+48)
     &        //CHAR(MOD(TIME_STEP,10)+48)
           END DO
      END IF

      IF (.NOT.USE_MPI) THEN
        OPEN(UNIT=10,FILE=FNAME,STATUS='NEW',FORM="UNFORMATTED")
        WRITE(10) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP

      IF (NUM_PER_DIR.EQ.3) THEN

         WRITE(10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     &             (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     &             (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY)

         DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN"
     &       ,FORM="UNFORMATTED")
          WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
          WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=0,TNKY)
          CLOSE(11)
         END DO
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        WRITE(10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY)
        DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN"
     &       ,FORM="UNFORMATTED")
          WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
          WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=1,NY)
          CLOSE(11)
        END DO
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        WRITE(10) (((CU1(I,K,J),I=0,NKX),K=1,NZ),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX),K=1,NZ),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX),K=2,NZ),J=1,NY)
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        WRITE(10) (((U1(I,K,J),I=2,NX),K=1,NZ),J=1,NY),
     *            (((U2(I,K,J),I=1,NX),K=1,NZ),J=2,NY),
     *            (((U3(I,K,J),I=1,NX),K=2,NZ),J=1,NY)
      END IF
      CLOSE(10)
      CLOSE(11)

      ELSE

        CALL MPI_IO_WRITE(FNAME,FNAME_TH,FNAME_P)
        IF (RANK.eq.0) then
           write(*,*) 'Saved flow to ',FNAME
           DO n=1,N_TH
             write(*,*) 'Saved flow to ',FNAME_TH(n)
           END DO
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      END IF

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FILTER(n)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      integer n

      IF (NUM_PER_DIR.EQ.3) THEN
        CALL FILTER_PER(n)
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        CALL FILTER_CHAN(n)
      END IF

      RETURN
      END

