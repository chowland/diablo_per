      program split_restart
! This rogram reads in a serial restart file and splits it into mupliple planes to be used by the MPI version
      include 'header'
      CHARACTER*35 FNAME
      CHARACTER*35 FNAME_TH(N_TH)
      integer NPROCS,i,j,k,n,np
   


      write(*,*) 'Enter the name of restart file: '
      read(*,*) FNAME
      do n=1,N_TH
        write(*,*) 'Enter the name of restart file for N=',n
        read(*,*) FNAME_TH(n)
      end do

      OPEN(UNIT=10,FILE=FNAME,STATUS="OLD",FORM="UNFORMATTED")
      READ (10) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP

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

        write(*,*) 'Enter the number of processors: '
        read(*,*) NPROCS

          do np=0,NPROCS-1

! First, get the character version of N
          IF (NPROCS.le.10) THEN
          MPI_IO_NUM=CHAR(MOD(np+1,10)+48)
        ELSE IF (NPROCS.le.100) THEN
          MPI_IO_NUM=CHAR(MOD(np+1,100)/10+48)
     &             //CHAR(MOD(np+1,10)+48)
        ELSE IF (NPROCS.le.1000) THEN
          MPI_IO_NUM=CHAR(MOD(np+1,1000)/100+48)
     &             //CHAR(MOD(np+1,100)/10+48)
     &             //CHAR(MOD(np+1,10)+48)
        ELSE IF (NPROCS.le.10000) THEN
          MPI_IO_NUM=CHAR(MOD(np+1,10000)/1000+48)
     &             //CHAR(MOD(np+1,1000)/100+48)
     &             //CHAR(MOD(np+1,100)/10+48)
     &             //CHAR(MOD(np+1,10)+48)
        ELSE
           WRITE(6,*) 'ERROR, NPROCS>10,000, Unsupported problem size'
        END IF


           FNAME='diablo_'//MPI_IO_NUM//'.start'
           DO N=1,N_TH
           FNAME_TH(N)='diablo_th'
     &        //CHAR(MOD(N,100)/10+48)
     &        //CHAR(MOD(N,10)+48) //'_'//MPI_IO_NUM//'.start'
           END DO

       WRITE(10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY)
        DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH(N),STATUS="UNKNOWN"
     &       ,FORM="UNFORMATTED")
          WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
          WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=0,NY)
          CLOSE(11)
        END DO

        end do
 
        stop
        end


  

