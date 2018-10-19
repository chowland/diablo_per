

!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header'

      INTEGER DIMS(2), PERDIM(2)
      INTEGER COMM_CART
! This subroutine initializes all mpi variables

      CALL MPI_INIT(IERROR)

      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERROR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,RANK,IERROR)

      NP_S=NPROCS**0.5
      RANKZ=mod(RANK,NP_S)
      RANKY=(RANK-RANKZ)/NP_S

      DIMS=(/NP_S,NP_S/)
      PERDIM=(/1,1/)

      call MPI_CART_CREATE(MPI_COMM_WORLD,2,DIMS,PERDIM,.FALSE.
     &                      ,COMM_CART,IERROR)
! In PERDIM, we put the information for the remain_dims
      PERDIM=(/1,0/)
      call MPI_CART_SUB(COMM_CART,PERDIM,MPI_COMM_Y,IERROR)
      PERDIM=(/0,1/)
      call MPI_CART_SUB(COMM_CART,PERDIM,MPI_COMM_Z,IERROR)

      call MPI_COMM_RANK(MPI_COMM_Y,RANKY,IERROR)
      call MPI_COMM_RANK(MPI_COMM_Z,RANKZ,IERROR)

      if (RANK.EQ.0) write(*,*) 'MPI Initialized. NPROCS: ',NPROCS
      call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      write(*,*) 'RANK,RANKY,RANKZ: ', RANK,RANKY,RANKZ

C Set a string to determine which input/output files to use
C When MPI is used, each process will read/write to files with the
C number of their rank (+1) appended to the end of the file.
C The string MPI_IO_NUM will be used to define the RANK+1 for each process
        IF (NPROCS.le.10) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,10)+48)
        ELSE IF (NPROCS.le.100) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,100)/10+48)
     &             //CHAR(MOD(RANK+1,10)+48)
        ELSE IF (NPROCS.le.1000) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,1000)/100+48)
     &             //CHAR(MOD(RANK+1,100)/10+48)
     &             //CHAR(MOD(RANK+1,10)+48)
        ELSE IF (NPROCS.le.10000) THEN
          MPI_IO_NUM=CHAR(MOD(RANK+1,10000)/1000+48)
     &             //CHAR(MOD(RANK+1,1000)/100+48)
     &             //CHAR(MOD(RANK+1,100)/10+48)
     &             //CHAR(MOD(RANK+1,10)+48)
        ELSE
           WRITE(6,*) 'ERROR, NPROCS>10,000, Unsupported problem size'
        END IF
        MPI_NUM=RANK

        NXM_S=NX_S-1
        NYM_S=NY_S-1
        NZM_S=NZ_S-1
        call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      SUBROUTINE TRANSPOSE_MPI_Z_TO_X(CA_Z,CA_X)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header'
! This subroutine initializes all mpi variables


      complex*16 CA_Z(0:NX_S/2,0:NZ+1,0:NY_S)
      complex*16 CA_X(0:NX/2,0:NZ_S,0:NY_S)

      complex*16 SENDBUF(0:NKX_S,0:NZ_S,0:NY_S)
      complex*16 RECVBUF(0:NKX_S,0:NZ_S,0:NY_S)

      integer DIAG,I,J,K
      integer block
      integer sendcount,recvcount
      integer shift,dir,dir_min
      integer I1,I2,J1,J2,K1,K2
      integer ISHIFT,JSHIFT,KSHIFT

! Zero Output array
      CA_X=(0.d0,0.d0)
      SENDBUF=(0.d0,0.d0)
      RECVBUF=(0.d0,0.d0)

      DIAG=MOD(RANK,NP_S)

      sendcount=(NKX_S+1)*(NZ_S+1)*(NY_S+1)
      recvcount=sendcount

      do shift=0,NP_S/2

      if ((shift*2.le.(NP_S-1)).and.(shift.ne.0)) then
        dir_min=-1
      else
        dir_min=1
      end if

      do dir=dir_min,1,2

      IF (SHIFT.eq.0) THEN
        BLOCK=DIAG
      ELSE
        BLOCK=DIAG + shift*dir*(1-2*MOD(DIAG/shift,2))
      END IF

      IF (BLOCK.GT.NP_S-1) THEN
        BLOCK=MOD(BLOCK,NP_S)
      ELSE IF (BLOCK.LT.0) THEN
        BLOCK=BLOCK+NP_S
      END IF

! Get a block from the source array to be sent to another process
        SENDBUF=0.d0
        K1=(NZ_S+1)*BLOCK
        K2=NZ+1-K1
        DO J=0,NY_S
          DO K=0,MIN(K2,NZ_S)
            DO I=0,NKX_S
              SENDBUF(I,K,J)=CA_Z(I,K+K1,J)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself
          RECVBUF=(0.d0,0.d0)
          CALL MPI_SENDRECV(
     &        sendbuf(0,0,0),sendcount,MPI_DOUBLE_COMPLEX
     &               ,block+(RANK/NP_S)*NP_S,12
     &        ,recvbuf(0,0,0),recvcount,MPI_DOUBLE_COMPLEX
     &               ,block+(RANK/NP_S)*NP_S,12
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        I1=(NKX_S+1)*BLOCK
        I2=NKX_S+I1
        DO J=0,NY_S
          DO K=0,NZ_S
            DO I=I1,MIN(I2,NX/2)
              CA_X(I,K,J)=RECVBUF(I-I1,K,J)
            END DO
          END DO
        END DO

        ELSE
! Here block=diag, we just need to copy locally from A->B
          IF (RANK.eq.0) THEN
! If we are on RANK=0, just copy directly from input to output
            DO J=0,NY_S
              DO K=0,NZ_S
                DO I=0,NKX_S
                  CA_X(I,K,J)=CA_Z(I,K,J)
                END DO
              END DO
            END DO
          ELSE
! Here, RANK.ne.0, so we need to shift the input and output
            K1=(NZ_S+1)*BLOCK
            K2=NZ+1-K1
            I1=(NKX_S+1)*BLOCK
            I2=I1+NX_S/2
            DO J=0,NY_S
              DO K=0,MIN(K2,NZ_S)
                DO I=I1,MIN(I2,NX/2)
                  CA_X(I,K,J)=CA_Z(I-I1,K+K1,J)
                END DO
              END DO
            END DO
          END IF

        END IF

      END DO

      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      SUBROUTINE TRANSPOSE_MPI_X_TO_Z(CA_X,CA_Z)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header'
! This subroutine initializes all mpi variables

      COMPLEX*16 CA_Z(0:NKX_S,0:NZ+1,0:NY_S)
      COMPLEX*16 CA_X(0:NX/2,0:NZ_S,0:NY_S)

      COMPLEX*16 SENDBUF(0:NKX_S,0:NZ_S,0:NY_S)
      COMPLEX*16 RECVBUF(0:NKX_S,0:NZ_S,0:NY_S)

      integer DIAG,I,J,K
      integer block
      integer shift,dir,dir_min
      integer sendcount,recvcount
      integer I1,I2,J1,J2,K1,K2
      integer ISHIFT,JSHIFT,KSHIFT

! Zero Output array
      CA_Z=(0.d0,0.d0)
      SENDBUF=(0.d0,0.d0)
      RECVBUF=(0.d0,0.d0)

      DIAG=MOD(RANK,NP_S)

      sendcount=(NKX_S+1)*(NZ_S+1)*(NY_S+1)
      recvcount=sendcount

      do shift=0,NP_S/2

      if ((shift*2.le.(NP_S-1)).and.(shift.ne.0)) then
        dir_min=-1
      else
        dir_min=1
      end if

      do dir=dir_min,1,2

      IF (SHIFT.eq.0) THEN
        BLOCK=DIAG
      ELSE
        BLOCK=DIAG + shift*dir*(1-2*MOD(DIAG/shift,2))
      END IF

      IF (BLOCK.GT.NP_S-1) THEN
        BLOCK=MOD(BLOCK,NP_S)
      ELSE IF (BLOCK.LT.0) THEN
        BLOCK=BLOCK+NP_S
      END IF

! Get a block from the source array to be sent to another process
        SENDBUF=0.d0
        I1=(NKX_S+1)*BLOCK
        I2=NX/2-I1
        DO J=0,NY_S
          DO K=0,NZ_S
            DO I=0,MIN(I2,NKX_S)
              SENDBUF(I,K,J)=CA_X(I+I1,K,J)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself
          RECVBUF=(0.d0,0.d0)
          CALL MPI_SENDRECV(
     &        sendbuf(0,0,0),sendcount,MPI_DOUBLE_COMPLEX
     &               ,block+(RANK/NP_S)*NP_S,13
     &        ,recvbuf(0,0,0),recvcount,MPI_DOUBLE_COMPLEX
     &               ,block+(RANK/NP_S)*NP_S,13
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        K1=(NZ_S+1)*BLOCK
        K2=NZ_S+K1
        DO J=0,NY_S
          DO K=K1,MIN(K2,NZ+1)
            DO I=0,NKX_S
! Make sure that we don't write beyond the end of the array
              CA_Z(I,K,J)=RECVBUF(I,K-K1,J)
            END DO
          END DO
        END DO

      ELSE
! Here block=diag, we just need to copy locally from A->B
        IF (RANK.eq.0) THEN
          DO J=0,NY_S
            DO K=0,NZ_S
              DO I=0,NKX_S
                CA_Z(I,K,J)=CA_X(I,K,J)
              END DO
            END DO
          END DO
        ELSE
          I1=(NKX_S+1)*BLOCK
          K1=(NZ_S+1)*BLOCK
          K2=NZ_S+K1
          DO J=0,NY_S
            DO K=K1,MIN(K2,NZ+1)
              DO I=0,NKX_S
                CA_Z(I,K,J)=CA_X(I+I1,K-K1,J)
              END DO
            END DO
          END DO
        END IF

        END IF

      END DO

      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      SUBROUTINE TRANSPOSE_MPI_X_TO_Y(A_X,A_Y,NX,NY,NZ,NX_S,NY_S,NZ_S)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables

      integer NX,NX_S,NY,NY_S,NZ,NZ_S
      complex*16 A_X(0:NX,0:NZ_S,0:NY_S)
      complex*16 A_Y(0:NX_S,0:NZ_S,0:NY)
      complex*16 B(0:NX_S,0:NZ_S,0:NY)

      complex*16 SENDBUF(0:NX_S,0:NZ_S,0:NY_S)
      complex*16 RECVBUF(0:NX_S,0:NZ_S,0:NY_S)

      integer DIAG,I,J,K
      integer block
      integer sendcount,recvcount

      DIAG=MOD(RANK,NP_S)

      sendcount=(NX_S+1)*(NZ_S+1)*(NY_S+1)
      recvcount=sendcount

      DO block=0,NP_S-1

! Get a block from the source array to be sent to another process
        SENDBUF=0.d0
        DO J=0,NY_S
          DO K=0,NZ_S
            DO I=0,NX_S
              SENDBUF(I,K,J)=A_X(I+(NX_S+1)*BLOCK,K,J)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself
          RECVBUF=(0.d0,0.d0)
          CALL MPI_SENDRECV(
     &        sendbuf,sendcount,MPI_DOUBLE_COMPLEX
     &               ,RANK+block-diag,RANK+block
     &        ,recvbuf,recvcount,MPI_DOUBLE_COMPLEX
     &               ,RANK+block-diag,RANK+block
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        DO I=0,NX_S
          DO K=0,NZ_S
            DO J=0,NY_S
              IF (J+(NY_S+1)*BLOCK.le.NY) THEN
! Make sure that we don't write beyond the end of the array
                B(I,K,J+(NY_S+1)*BLOCK)=RECVBUF(I,K,J)
              END IF
            END DO
          END DO
        END DO

        ELSE
! Here block=diag, we just need to copy locally from A->B



          DO I=0,NX_S
            DO K=0,NZ_S
              DO J=0,NY_S
                sendbuf(I,K,J)=A_X(I+(NX_S+1)*BLOCK,K,J)
              END DO
            END DO
          END DO

          DO I=0,NX_S
            DO K=0,NZ_S
              DO J=0,NY_S
              IF (J+(NY_S+1)*BLOCK.le.NY) THEN
                B(I,K,J+(NY_S+1)*BLOCK)=sendbuf(I,K,J)
              END IF
              END DO
            END DO
          END DO
        END IF

      END DO

! Finally, copy the temporary storage array to the final array
      DO I=0,NX_S
        DO K=0,NZ_S
          DO J=0,NY
            A_Y(I,K,J)=B(I,K,J)
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      SUBROUTINE TRANSPOSE_MPI_Y_TO_X(A_Y,A_X,NX,NY,NZ,NX_S,NY_S,NZ_S)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables

      integer NX,NX_S,NY,NY_S,NZ,NZ_S
      complex*16 A_X(0:NX,0:NZ_S,0:NY_S)
      complex*16 A_Y(0:NX_S,0:NZ_S,0:NY)
      complex*16 B(0:NX,0:NZ_S,0:NY_S)

      complex*16 SENDBUF(0:NX_S,0:NZ_S,0:NY_S)
      complex*16 RECVBUF(0:NX_S,0:NZ_S,0:NY_S)

      integer DIAG,I,J,K
      integer block
      integer sendcount,recvcount

      integer I1,I2,J1,J2,K1,K2

      DIAG=MOD(RANK,NP_S)

      sendcount=(NX_S+1)*(NZ_S+1)*(NY_S+1)
      recvcount=sendcount

      DO block=0,NP_S-1

! Get a block from the source array to be sent to another process
        SENDBUF=0.d0
        DO I=0,NX_S
          DO K=0,NZ_S
            DO J=0,NY_S
              SENDBUF(I,K,J)=A_Y(I,K,J+(NY_S+1)*BLOCK)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
          RECVBUF=(0.d0,0.d0)
! Communicate only if we aren't talking to ourself
          CALL MPI_SENDRECV(
     &        sendbuf,sendcount,MPI_DOUBLE_COMPLEX
     &               ,RANK+block-diag,RANK+block
     &        ,recvbuf,recvcount,MPI_DOUBLE_COMPLEX
     &               ,RANK+block-diag,RANK+block
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        I2=NX-BLOCK*(NX_S+1)
        DO I=0,MIN(I2,NX_S)
          DO K=0,NZ_S
            DO J=0,NY_S
! Make sure that we don't write beyond the end of the array
              B(I+(NX_S+1)*BLOCK,K,J)=RECVBUF(I,K,J)
            END DO
          END DO
        END DO

        ELSE
! Here block=diag, we just need to copy locally from A->B

          DO I=0,NX_S
            DO K=0,NZ_S
              DO J=0,NY_S
                sendbuf(I,K,J)=A_Y(I,K,J+(NY_S+1)*BLOCK)
              END DO
            END DO
          END DO
          DO I=0,NX_S
            DO K=0,NZ_S
              DO J=0,NY_S
              IF (I+(NX_S+1)*BLOCK.le.NX) THEN
                B(I+(NX_S+1)*BLOCK,K,J)=sendbuf(I,K,J)
              END IF
              END DO
            END DO
          END DO
        END IF

      END DO

! Finally, copy the temporary storage array to the final array
      DO I=0,NX
        DO K=0,NZ_S
          DO J=0,NY_S
            A_X(I,K,J)=B(I,K,J)
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      SUBROUTINE TRANSPOSE_MPI_Y_TO_Z(CA_Y,CA_Z)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header'
! This subroutine initializes all mpi variables


      complex*16 CA_Z(0:NX_S/2,0:NZ+1,0:NY_S)
      complex*16 CA_Y(0:NX_S/2,0:NZ_S,0:NY+1)

      complex*16 SENDBUF(0:NKX_S,0:TNKZ_S,0:NY_S)
      complex*16 RECVBUF(0:NKX_S,0:TNKZ_S,0:NY_S)

      integer DIAG,I,J,K
      integer block
      integer sendcount,recvcount
      integer shift,dir,dir_min
      integer I1,I2,K1,K2,J1,J2
      integer ISHIFT,JSHIFT,KSHIFT

! Zero Output array
      CA_Z=(0.d0,0.d0)
      SENDBUF=(0.d0,0.d0)
      RECVBUF=(0.d0,0.d0)

      DIAG=RANK/NP_S

      sendcount=(NKX_S+1)*(TNKZ_S+1)*(NY_S+1)
      recvcount=sendcount

      do shift=0,NP_S/2


      if ((shift*2.le.(NP_S-1)).and.(shift.ne.0)) then
        dir_min=-1
      else
        dir_min=1
      end if

      do dir=dir_min,1,2

      IF (SHIFT.eq.0) THEN
        BLOCK=DIAG
      ELSE
        BLOCK=DIAG + shift*dir*(1-2*MOD(DIAG/shift,2))
      END IF

      IF (BLOCK.GT.NP_S-1) THEN
        BLOCK=MOD(BLOCK,NP_S)
      ELSE IF (BLOCK.LT.0) THEN
        BLOCK=BLOCK+NP_S
      END IF

! Get a block from the source array to be sent to another process
        SENDBUF=0.d0
        J1=(NY_S+1)*BLOCK
        J2=NY+1-J1
        DO J=0,MIN(J2,NY_S)
          DO K=0,TNKZ_S
            DO I=0,NKX_S
              SENDBUF(I,K,J)=CA_Y(I,K,J+J1)
            END DO
          END DO
        END DO


        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself
      RECVBUF=(0.d0,0.d0)
          CALL MPI_SENDRECV(
     &        sendbuf,sendcount,MPI_DOUBLE_COMPLEX
     &               ,MOD(RANK,NP_S)+block*NP_S
     &               ,RANK+MOD(RANK,NP_S)+block*NP_S
     &        ,recvbuf,recvcount,MPI_DOUBLE_COMPLEX
     &               ,MOD(RANK,NP_S)+block*NP_S
     &               ,RANK+MOD(RANk,NP_S)+block*NP_S
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        K1=(TNKZ_S+1)*BLOCK
        K2=TNKZ_S+K1
        DO J=0,NY_S
          DO K=K1,MIN(K2,NZ+1)
            DO I=0,NKX_S
              CA_Z(I,K,J)=RECVBUF(I,K-K1,J)
            END DO
          END DO
        END DO

        ELSE
! Here block=diag, we just need to copy locally from A->B
          IF (RANK.eq.0) THEN
          DO J=0,NY_S
            DO K=0,TNKZ_S
              DO I=0,NKX_S
                CA_Z(I,K,J)=CA_Y(I,K,J)
              END DO
            END DO
          END DO
          ELSE
          J1=(NY_S+1)*BLOCK
          J2=NY+1-J1
          K1=(TNKZ_S+1)*BLOCK
          K2=K1+TNKZ_S
          DO J=0,MIN(J2,NY_S)
            DO K=K1,MIN(K2,NZ+1)
              DO I=0,NKX_S
                CA_Z(I,K,J)=CA_Y(I,K-K1,J+J1)
              END DO
            END DO
          END DO
          END IF

        END IF


      END DO

      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      SUBROUTINE TRANSPOSE_MPI_Z_TO_Y(CA_Z,CA_Y)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header'
! This subroutine initializes all mpi variables

      complex*16 CA_Z(0:NKX_S,0:NZ+1,0:NY_S)
      complex*16 CA_Y(0:NKX_S,0:TNKZ_S,0:NY+1)

      complex*16 SENDBUF(0:NKX_S,0:TNKZ_S,0:NY_S)
      complex*16 RECVBUF(0:NKX_S,0:TNKZ_S,0:NY_S)

      integer DIAG,I,J,K
      integer block
      integer sendcount,recvcount
      integer shift,dir,dir_min
      integer I1,I2,J1,J2,K1,K2
      integer ISHIFT,JSHIFT,KSHIFT

! Zero Output array, Working arrays
      CA_Y=(0.d0,0.d0)
      SENDBUF=(0.d0,0.d0)
      RECVBUF=(0.d0,0.d0)

      DIAG=RANK/NP_S

      sendcount=(NKX_S+1)*(TNKZ_S+1)*(NY_S+1)
      recvcount=sendcount

      do shift=0,NP_S/2

      if ((shift*2.le.(NP_S-1)).and.(shift.ne.0)) then
        dir_min=-1
      else
        dir_min=1
      end if

      do dir=dir_min,1,2

      IF (SHIFT.eq.0) THEN
         BLOCK=DIAG
      ELSE
        BLOCK=DIAG + shift*dir*(1-2*MOD(DIAG/shift,2))
      END IF

      IF (BLOCK.GT.NP_S-1) THEN
        BLOCK=MOD(BLOCK,NP_S)
      ELSE IF (BLOCK.LT.0) THEN
        BLOCK=BLOCK+NP_S
      END IF

! Get a block from the source array to be sent to another process
        SENDBUF=0.d0
        K1=(TNKZ_S+1)*BLOCK
        K2=NZ+1-K1
        DO J=0,NY_S
          DO K=0,MIN(K2,TNKZ_S)
            DO I=0,NKX_S
              SENDBUF(I,K,J)=CA_Z(I,K+K1,J)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
          RECVBUF=(0.d0,0.d0)
! Communicate only if we aren't talking to ourself
          CALL MPI_SENDRECV(
     &        sendbuf,sendcount,MPI_DOUBLE_COMPLEX
     &               ,MOD(RANK,NP_S)+block*NP_S,14
     &        ,recvbuf,recvcount,MPI_DOUBLE_COMPLEX
     &               ,MOD(RANK,NP_S)+block*NP_S,14
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        J1=(NY_S+1)*BLOCK
        J2=NY_S+J1
        DO J=J1,MIN(J2,NY+1)
          DO K=0,TNKZ_S
            DO I=0,NKX_S
! Make sure that we don't write beyond the end of the array
               CA_Y(I,K,J)=RECVBUF(I,K,J-J1)
            END DO
          END DO
        END DO

        ELSE
! Here block=diag, we just need to copy locally from A->B
          IF (RANK.eq.0) THEN
          DO J=0,NY_S
              DO K=0,TNKZ_S
                DO I=0,NKX_S
                  CA_Y(I,K,J)=CA_Z(I,K,J)
                END DO
              END DO
          END DO
          ELSE
          K1=(TNKZ_S+1)*BLOCK
          K2=NZ+1-K1
          J1=(NY_S+1)*BLOCK
          J2=NY_S+J1
          DO J=J1,MIN(J2,NY+1)
              DO K=0,MIN(K2,TNKZ_S)
                DO I=0,NKX_S
                  CA_Y(I,K,J)=CA_Z(I,K+K1,J-J1)
                END DO
              END DO
          END DO
          END IF
        END IF

      END DO

      END DO


      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      subroutine courant_mpi
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
! This subroutine sets the timestep based on the specified CFL number
! The subroutine should be called with the velocity in physical space

      include 'header'

      real*8 vel
      real*8 dt_local,dt_test
      real*8 dt_x,dt_y,dt_z
      real*8 dt
      integer i,j,k
      integer imin,jmin,kmin

! Set the initial dt to some arbitrary large number
! For chemotaxis with flow
      dt_local=7.d-3
! For no flow: (rescaled for 128^3)
!      dt_local=7.d-2*25.d0
! Set the maximum dt based on the buoyancy period
!      dt_local=0.005d0*2.d0*PI/sqrt(RI_TAU(1))

      dt_test=0.d0

      do j=0,NY_S
        do k=0,NZ_S
          do i=0,NXM
            IF (U1(i,k,j).ne.0.d0) then
              dt_x=abs(cfl*DX_TH(i)/abs(U1(i,k,j)))
            ELSE
              dt_x=999.d0
            END IF
            IF (U2(i,k,j).ne.0.d0) then
              dt_y=abs(cfl*DY_TH(j)/abs(U2(i,k,j)))
            ELSE
              dt_y=999.d0
            END IF
! Note, thermal wind advection included
!            dt_z=abs(cfl*dz(k)/abs(
!     *     U3(i,k,j)+(RI_TAU(N)/I_RO_TAU)*DRHODX*GYF(J)
!     &      ))
            IF (U3(i,k,j).ne.0.d0) THEN
              dt_z=abs(cfl*DZ_TH(k)/abs(U3(i,k,j)))
            ELSE
              dt_z=999.d0
            END IF
            dt_local=min(dt_local,dt_x,dt_y,dt_z)
          end do
        end do
      end do

! Now we have the minimum DELTA_T on each process
! share the information to obtain the global minimum

      CALL MPI_ALLREDUCE(dt_local,dt,1
     &         ,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)

      if (dt.le.0) then
        write(*,*) 'Error: dt<=0 in courant'
! Set DELTA_T to some small default value
        DELTA_T=0.0001d0
      else if (dt.ge.1000.) then
        write(*,*) 'WARNING: DELTA_T > 1000, value capped at 1000'
        DELTA_T=1000.d0
        H_BAR(1)=DELTA_T*(8.0/15.0)
        H_BAR(2)=DELTA_T*(2.0/15.0)
        H_BAR(3)=DELTA_T*(5.0/15.0)
      else
        DELTA_T=dt
        H_BAR(1)=DELTA_T*(8.0/15.0)
        H_BAR(2)=DELTA_T*(2.0/15.0)
        H_BAR(3)=DELTA_T*(5.0/15.0)
      end if

      return
      end


C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      SUBROUTINE END_RUN_MPI(FLAG)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header'

      LOGICAL FLAG

      IF (RANK.EQ.0) THEN
        CALL END_RUN(FLAG)
      END IF
      CALL MPI_BCAST(FLAG,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERROR)

      END