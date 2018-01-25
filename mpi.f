

!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header'

      INTEGER DIMS(2), PERDIM(2)
! This subroutine initializes all mpi variables

      CALL MPI_INIT(IERROR)
      write(*,*) 'MPI_COMM_WORLD: ',MPI_COMM_WORLD

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


      SUBROUTINE TRANSPOSE_MPI_Z_TO_X(CA_Z,CA_X)
      INCLUDE 'header'
! This subroutine initializes all mpi variables


      complex*16 CA_Z(0:NX_S/2,0:NZ+1,0:NY_S)
      complex*16 CA_X(0:NX/2,0:NZ_S,0:NY_S)

      complex*16 SENDBUF(0:NKX_S,0:NZ_S,0:NY_S)
      complex*16 RECVBUF(0:NKX_S,0:NZ_S,0:NY_S)

      integer DIAG
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

      SUBROUTINE TRANSPOSE_MPI_X_TO_Z(CA_X,CA_Z)
      INCLUDE 'header'
! This subroutine initializes all mpi variables

      COMPLEX*16 CA_Z(0:NKX_S,0:NZ+1,0:NY_S)
      COMPLEX*16 CA_X(0:NX/2,0:NZ_S,0:NY_S)

      COMPLEX*16 SENDBUF(0:NKX_S,0:NZ_S,0:NY_S)
      COMPLEX*16 RECVBUF(0:NKX_S,0:NZ_S,0:NY_S)

      integer DIAG
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

      SUBROUTINE TRANSPOSE_MPI_X_TO_Y(A_X,A_Y,NX,NY,NZ,NX_S,NY_S,NZ_S)
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables


      complex*16 A_X(0:NX,0:NZ_S,0:NY_S)
      complex*16 A_Y(0:NX_S,0:NZ_S,0:NY)
      complex*16 B(0:NX_S,0:NZ_S,0:NY)

      complex*16 SENDBUF(0:NX_S,0:NZ_S,0:NY_S)
      complex*16 RECVBUF(0:NX_S,0:NZ_S,0:NY_S)

      integer DIAG
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

      SUBROUTINE TRANSPOSE_MPI_Y_TO_X(A_Y,A_X,NX,NY,NZ,NX_S,NY_S,NZ_S)
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables


      complex*16 A_X(0:NX,0:NZ_S,0:NY_S)
      complex*16 A_Y(0:NX_S,0:NZ_S,0:NY)
      complex*16 B(0:NX,0:NZ_S,0:NY_S)

      complex*16 SENDBUF(0:NX_S,0:NZ_S,0:NY_S)
      complex*16 RECVBUF(0:NX_S,0:NZ_S,0:NY_S)

      integer DIAG
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

      SUBROUTINE TRANSPOSE_MPI_Y_TO_Z(CA_Y,CA_Z)
      INCLUDE 'header'
! This subroutine initializes all mpi variables


      complex*16 CA_Z(0:NX_S/2,0:NZ+1,0:NY_S)
      complex*16 CA_Y(0:NX_S/2,0:NZ_S,0:NY+1)

      complex*16 SENDBUF(0:NKX_S,0:TNKZ_S,0:NY_S)
      complex*16 RECVBUF(0:NKX_S,0:TNKZ_S,0:NY_S)

      integer DIAG
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

      SUBROUTINE TRANSPOSE_MPI_Z_TO_Y(CA_Z,CA_Y)
      INCLUDE 'header'
! This subroutine initializes all mpi variables

      complex*16 CA_Z(0:NKX_S,0:NZ+1,0:NY_S)
      complex*16 CA_Y(0:NKX_S,0:TNKZ_S,0:NY+1)

      complex*16 SENDBUF(0:NKX_S,0:TNKZ_S,0:NY_S)
      complex*16 RECVBUF(0:NKX_S,0:TNKZ_S,0:NY_S)

      integer DIAG
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





      subroutine courant_mpi
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


!      IF (RANK.eq.0) THEN
!       CALL MPI_SEND(dt,1,MPI_DOUBLE_PRECISION,RANK+1,1
!     &     ,MPI_COMM_WORLD
!     &     ,IERROR)
!
!      ELSE IF (RANK.lt.NPROCS-1) THEN
!       CALL MPI_RECV(dt_test,1,MPI_DOUBLE_PRECISION,RANK-1,1
!     &     ,MPI_COMM_WORLD
!     &     ,STATUS,IERROR)
!       dt=min(dt,dt_test)
!       CALL MPI_SEND(dt,1,MPI_DOUBLE_PRECISION,RANK+1,1
!     &     ,MPI_COMM_WORLD
!     &     ,IERROR)
!       ELSE IF (RANK.eq.NPROCS-1) THEN
!        CALL MPI_RECV(dt_test,1,MPI_DOUBLE_PRECISION,RANK-1,1
!     &     ,MPI_COMM_WORLD
!     &     ,STATUS,IERROR)
!        dt=min(dt,dt_test)
!      END IF
! Now, the process with RANK=NPROCS-1 has the minimum dt
! share this with the other processes
!      IF (RANK.eq.NPROCS-1) THEN
!       CALL MPI_SEND(dt,1,MPI_DOUBLE_PRECISION,RANK-1,1
!     &     ,MPI_COMM_WORLD
!     &     ,IERROR)
!      ELSE IF (RANK.gt.0) THEN
!        CALL MPI_RECV(dt,1,MPI_DOUBLE_PRECISION,RANK+1,1
!     &     ,MPI_COMM_WORLD
!     &     ,STATUS,IERROR)
!        CALL MPI_SEND(dt,1,MPI_DOUBLE_PRECISION,RANK-1,1
!     &     ,MPI_COMM_WORLD
!     &     ,IERROR)
!      ELSE
!        CALL MPI_RECV(dt,1,MPI_DOUBLE_PRECISION,RANK+1,1
!     &     ,MPI_COMM_WORLD
!     &     ,STATUS,IERROR)
!      END IF


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




      SUBROUTINE GHOST_CHAN_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep

      include 'header'

      integer i,j,k,N

! Define the arrays that will be used for data packing.  This makes the
! communication between processes more efficient by only requiring one
! send and recieve.
! The communication will be done in Fourier space, so these arrays should
! be complex arrays to match the velocity
! The size of the buffer arrys is 0:NKX,0:TNKZ,# of variables
      COMPLEX*16 OCPACK(0:NX/3,0:2*NZ/3,4+N_TH)
      COMPLEX*16 ICPACK(0:NX/3,0:2*NZ/3,4+N_TH)

! If we are using more than one processor, then we need to pass data
      IF (NPROCS.gt.1) THEN

! First, Pass data up the chain to higher ranked processes

      IF (RANK.eq.0) THEN
! If we are the lowest ranked process, then we don't need to recieve
! data at the lower ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)


! End if RANK=0
      ELSE IF (RANK.lt.NPROCS-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,NY)
            OCPACK(I,K,2)=CU2(I,K,NY)
            OCPACK(I,K,3)=CU3(I,K,NY)
            OCPACK(I,K,4)=CP(I,K,NY)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,NY,N)
            END DO
          END DO
        END DO
! Use MPI_SENDRECV since we need to recieve and send data
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved

        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,1)=ICPACK(I,K,1)
            CU2(I,K,1)=ICPACK(I,K,2)
            CU3(I,K,1)=ICPACK(I,K,3)
            CP(I,K,1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF

! AT this point we have passed data up the chain

! Now, we need to pass data down the chain of processes
      IF (RANK.eq.NPROCS-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO
! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
      ELSE IF (RANK.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us
        DO K=0,TNKZ
          DO I=0,NKX
            OCPACK(I,K,1)=CU1(I,K,2)
            OCPACK(I,K,2)=CU2(I,K,2)
            OCPACK(I,K,3)=CU3(I,K,2)
            OCPACK(I,K,4)=CP(I,K,2)
            DO N=1,N_TH
              OCPACK(I,K,4+N)=CTH(I,K,2,N)
            END DO
          END DO
        END DO

        CALL MPI_SEND(OCPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,(4+N_TH)*(NKX+1)*(TNKZ+1)
     &               ,MPI_DOUBLE_COMPLEX
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,NY+1)=ICPACK(I,K,1)
            CU2(I,K,NY+1)=ICPACK(I,K,2)
            CU3(I,K,NY+1)=ICPACK(I,K,3)
            CP(I,K,NY+1)=ICPACK(I,K,4)
            DO N=1,N_TH
              CTH(I,K,NY+1,N)=ICPACK(I,K,4+N)
            END DO
          END DO
        END DO
      END IF

      END IF

      RETURN
      END


      SUBROUTINE GHOST_GRID_MPI
! This subroutine is part of the MPI package for the channel flow
! Diablo package.
! Here, we define a set of ghost cells on each process
! the ghost cells contain information from the neighboring nodes
! and allow us to compute finite differences over the local gridpoints.
! We need to update the contents of the ghost cells at the start of
! each Runge-Kutta substep

      include 'header'

      integer i,j,k,N

      real*8 OCPACK(3),ICPACK(3)

! First, Pass data up the chain to higher ranked processes

      IF (RANK.eq.0) THEN

! Set the lower ghost cells
        GYF(0)=2.d0*GYF(1)-GYF(2)
        GY(0)=2.d0*GY(1)-GY(2)

        OCPACK(1)=GYF(NY-1)
        OCPACK(2)=GYF(NY)
        OCPACK(3)=GY(NY)

        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

! End if RANK=0
      ELSE IF (RANK.lt.NPROCS-1) THEN
! Here, we are one of the middle processes and we need to pass data
! up and recieve data from below
        OCPACK(1)=GYF(NY-1)
        OCPACK(2)=GYF(NY)
        OCPACK(3)=GY(NY)

        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,1,MPI_COMM_WORLD,IERROR)

        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(0)=ICPACK(1)
        GYF(1)=ICPACK(2)
        GY(1)=ICPACK(3)

      ELSE
! Otherwise, we must be the uppermost process with RANK=NPROCS-1
! Set the top ghost cell
        GYF(NY+1)=2.*GYF(NY)-GYF(NYM)

! Here, we need to recieve data from below, but don't need to send data up
        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,1,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(0)=ICPACK(1)
        GYF(1)=ICPACK(2)
        GY(1)=ICPACK(3)

      END IF
! AT this point we have passed data up the chain

      IF (RANK.eq.NPROCS-1) THEN
! If we are the higest ranked process, then we don't need to recieve
! data at the upper ghost cells, these will be filled with boundary
! condition information

        OCPACK(1)=GYF(2)
        OCPACK(2)=GY(2)


! Now, we have packed the data into a compact array, pass the data up
        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
      ELSE IF (RANK.GT.0) THEN
! Here, we are one of the middle processes and we need to pass data
! down and recieve data from above us

        OCPACK(1)=GYF(2)
        OCPACK(2)=GY(2)
        CALL MPI_SEND(OCPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK-1,3,MPI_COMM_WORLD,IERROR)
        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Now, unpack the data that we have recieved
        GYF(NY+1)=ICPACK(1)
        GY(NY+1)=ICPACK(2)

      ELSE
! Here, we must be the lowest process (RANK=0) and we need to recieve
! data from above
        CALL MPI_RECV(ICPACK,3
     &               ,MPI_DOUBLE_PRECISION
     &               ,RANK+1,3,MPI_COMM_WORLD,STATUS,IERROR)
! Unpack the data that we have recieved
        GYF(NY+1)=ICPACK(1)
        GY(NY+1)=ICPACK(2)

      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE THOMAS_FORWARD_REAL_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This version solves for one full row, then passes the data
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX-1,0:NY+1), B(0:NX-1,0:NY+1), C(0:NX-1,0:NY+1)
     &         , G(0:NX-1,0:NY+1)

      REAL*8 OCPACK(0:NX-1,4),ICPACK(0:NX-1,4)


      IF (RANK.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4*NX,MPI_DOUBLE_PRECISION,RANK-1,12
     &               ,MPI_COMM_WORLD,status,ierror)
C Unpack the data
        DO I=0,NX-1
        A(I,1)=OCPACK(I,1)
        B(I,1)=OCPACK(I,2)
        C(I,1)=OCPACK(I,3)
        G(I,1)=OCPACK(I,4)
        END DO
      ELSE
! We are the lowest process, so solve at J=1
        J=1
        DO I=0,NX-1
          A(I,J)=-A(I,J)/B(I,J-1)
          B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
          G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
        END DO
       END IF

! For all processes, solve from 2..NY
        DO J=2,NY
          DO I=0,NX-1
            A(I,J)=-A(I,J)/B(I,J-1)
            B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
            G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
          END DO
        END DO

      IF (RANK.NE.NPROCS-1) THEN
        DO I=0,NX-1
        ICPACK(I,1)=A(I,NY)
        ICPACK(I,2)=B(I,NY)
        ICPACK(I,3)=C(I,NY)
        ICPACK(I,4)=G(I,NY)
        END DO
        CALL MPI_SEND(ICPACK,4*NX,MPI_DOUBLE_PRECISION,RANK+1,12
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_FORWARD_COMPLEX_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX,0:NY+1), B(0:NX,0:NY+1), C(0:NX,0:NY+1)
      COMPLEX*16  G(0:NX,0:NY+1)

      COMPLEX*16 OCPACK(0:NX,4),ICPACK(0:NX,4)


      IF (RANK.NE.0) THEN
C If we aren't the lowest process, then wait for data
        CALL MPI_RECV(OCPACK,4*(NX+1),MPI_DOUBLE_COMPLEX,RANK-1,13
     &               ,MPI_COMM_WORLD,status,ierror)
C Unpack the data
        DO I=0,NX
        A(I,1)=real(OCPACK(I,1))
        B(I,1)=real(OCPACK(I,2))
        C(I,1)=real(OCPACK(I,3))
        G(I,1)=OCPACK(I,4)
        END DO
       ELSE
! We are the lowest process, so solve at J=0
        J=1
        DO I=0,NX
        A(I,J)=-A(I,J)/B(I,J-1)
        B(I,J)=B(I,J)+A(I,J)*C(I,J-1)
        G(I,J)=G(I,J)+A(I,J)*G(I,J-1)
        END DO
       END IF

! For all processes, solve from 2..NY
        DO J=2,NY
          DO I=0,NX
            A(I,J)=-A(I,J)/B(I,J-1)
            B(I,J)=B(I,J)+A(I,J+1)*C(I,J-1)
            G(I,J)=G(I,J)+A(I,J+1)*G(I,J-1)
          END DO
        END DO

      IF (RANK.NE.NPROCS-1) THEN
        DO I=0,NX
        ICPACK(I,1)=A(I,NY)
        ICPACK(I,2)=B(I,NY)
        ICPACK(I,3)=C(I,NY)
        ICPACK(I,4)=G(I,NY)
        END DO
        CALL MPI_SEND(ICPACK,4*(NX+1),MPI_DOUBLE_COMPLEX,RANK+1,13
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      RETURN
      END




C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_REAL_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|------
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX-1,0:NY+1), B(0:NX-1,0:NY+1), C(0:NX-1,0:NY+1)
     &     , G(0:NX-1,0:NY+1)
      REAL*8 ICPACK(0:NX-1),OCPACK(0:NX-1)

      IF (RANK.NE.NPROCS-1) THEN
C If we aren't the highest process, then wait for data
        CALL MPI_RECV(OCPACK,NX,MPI_DOUBLE_PRECISION,RANK+1,10
     &               ,MPI_COMM_WORLD,status,ierror)
        DO I=0,NX-1
          G(I,NY+1)=OCPACK(I)
        END DO
      ELSE
C Else, if we are the highest process, compute the solution at j=NY
        DO I=0,NX-1
          G(I,NY)=G(I,NY)/B(I,NY)
        END DO
      END IF

C All processes solve from NY..1
      DO J=NY,1,-1
      DO I=0,NX-1
        G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
      END DO
      END DO

      IF (RANK.NE.0) THEN
        DO I=0,NX-1
          ICPACK(I)=G(I,2)
        END DO
        CALL MPI_SEND(ICPACK,NX,MPI_DOUBLE_PRECISION,RANK-1,10
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE THOMAS_BACKWARD_COMPLEX_MPI(A,B,C,G,NY,NX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C This subroutine performs the backward sweep of the Thomas algorithm
C Thomas algorithm solves Ax=b for tridiagonal A
C The RHS vector and solution are real
C Input lower, main, and upper diagonals, ld, md, ud, and rhs x
C Returns solution in x
C The indexing should be done by ROW, ie.
C [ b1  c1   0   0   0 ...
C [ a2  b2  c2   0   0 ...
C [  0  a3  b3   c3  0 ...

      INCLUDE 'header_mpi'

      INTEGER I,J,NX,NY
      REAL*8 A(0:NX,0:NY+1), B(0:NX,0:NY+1), C(0:NX,0:NY+1)
      COMPLEX*16 G(0:NX,0:NY+1)
      COMPLEX*16 ICPACK(0:NX),OCPACK(0:NX)

      IF (RANK.NE.NPROCS-1) THEN
C If we aren't the highest process, then wait for data
        CALL MPI_RECV(OCPACK,NX+1,MPI_DOUBLE_COMPLEX,RANK+1,11
     &               ,MPI_COMM_WORLD,status,ierror)
        DO I=0,NX
          G(I,NY+1)=OCPACK(I)
        END DO
      ELSE
C Else, if we are the highest process, compute the solution directly at j=NY
        DO I=0,NX
          G(I,NY)=G(I,NY)/B(I,NY)
        END DO
      END IF

C All processes solve from NY-1..0
      DO J=NY,1,-1
      DO I=0,NX
        G(I,J)=(G(I,J)-C(I,J)*G(I,J+1))/B(I,J)
      END DO
      END DO

      IF (RANK.NE.0) THEN
        DO I=0,NX
          ICPACK(I)=G(I,2)
        END DO
        CALL MPI_SEND(ICPACK,NX+1,MPI_DOUBLE_COMPLEX,RANK-1,11
     &               ,MPI_COMM_WORLD,ierror)
      END IF

      RETURN
      END






!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_CHAN_MPI
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
C Initialize any constants here
      INCLUDE 'header'

      INTEGER J, N

      write(*,*) '*******IN INIT_CHAN_MPI*********'


      PI=4.D0*ATAN(1.D0)

! Set the upper and lower bounds for timestepping
      IF (RANK.eq.0) THEN
        JEND=NY
        IF (U_BC_YMIN.EQ.0) THEN
          JSTART=2
        ELSE IF (U_BC_YMIN.EQ.1) THEN
          JSTART=1
        ELSE
          JSTART=2
        END IF
! Now, set the indexing for the scalar equations
        DO N=1,N_TH
          JEND_TH(N)=NY
          IF (TH_BC_YMIN(N).EQ.0) THEN
            JSTART_TH(N)=2
          ELSE IF (TH_BC_YMIN(N).EQ.1) THEN
            JSTART_TH(N)=1
          ELSE
            JSTART_TH(N)=2
          END IF
        END DO
      ELSE IF (RANK.eq.NPROCS-1) THEN
        JSTART=2
        IF (U_BC_YMAX.EQ.0) THEN
          JEND=NY-1
        ELSE IF (U_BC_YMAX.EQ.1) THEN
          JEND=NY
        ELSE
          JEND=NY-1
        END IF

! Set the upper and lower limits of timestepping of the scalar equations
        DO N=1,N_TH
        JSTART_TH(N)=2
        IF (TH_BC_YMAX(N).EQ.0) THEN
          JEND_TH(N)=NY-1
        ELSE IF (TH_BC_YMAX(N).EQ.1) THEN
          JEND_TH(N)=NY
        ELSE
          JEND_TH(N)=NY-1
        END IF
        END DO

      ELSE
! Here, we are on a middle process
        JSTART=2
        JEND=NY
        DO N=1,N_TH
          JSTART_TH(N)=2
          JEND_TH(N)=NY
        END DO
      END IF

      RETURN
      END


      SUBROUTINE APPLY_BC_TH_MPI(MATL,MATD,MATU,VEC,N)
! This subroutine applies the boundary conditions to the
! scalar fields prior to the implicit solve
      INCLUDE 'header'
      INTEGER N
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_TH_LOWER(MATL,MATD,MATU,VEC,N)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the upper plane, apply the boundary conditions
          CALL APPLY_BC_TH_UPPER(MATL,MATD,MATU,VEC,N)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_TH_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C,K,N)
! This subroutine applies the boundary conditions to the
! scalar fields prior to the implicit solve
      INCLUDE 'header'
      INTEGER N,K
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_TH_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C,K,N)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the upper plane, apply the boundary conditions
          CALL APPLY_BC_TH_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C,K,N)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_U2_MPI(MATL,MATD,MATU,VEC)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_2_LOWER(MATL,MATD,MATU,VEC)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions
          CALL APPLY_BC_2_UPPER(MATL,MATD,MATU,VEC)
        END IF
      RETURN
      END

      SUBROUTINE APPLY_BC_U2_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      integer K
      INCLUDE 'header'
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_2_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C,K)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions
          CALL APPLY_BC_2_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C,K)
        END IF
      RETURN
      END
      SUBROUTINE APPLY_BC_U1_MPI(MATL,MATD,MATU,VEC)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_1_LOWER(MATL,MATD,MATU,VEC)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions
          CALL APPLY_BC_1_UPPER(MATL,MATD,MATU,VEC)
        END IF
      RETURN
      END
      SUBROUTINE APPLY_BC_U1_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      integer k
      INCLUDE 'header'
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_1_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C,K)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions
          CALL APPLY_BC_1_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C,K)
        END IF
      RETURN
      END
      SUBROUTINE APPLY_BC_U3_MPI(MATL,MATD,MATU,VEC)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_3_LOWER(MATL,MATD,MATU,VEC)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions
          CALL APPLY_BC_3_UPPER(MATL,MATD,MATU,VEC)
        END IF
      RETURN
      END
      SUBROUTINE APPLY_BC_U3_MPI_C(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions to the
! velocity field prior to the implicit solve
! Note, MATL, MATD, etc. are dimensioned in header
      integer K
      INCLUDE 'header'
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
        IF (RANK.eq.0) THEN
! If we have the lowest plane, apply the boundary conditions
          CALL APPLY_BC_3_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C,K)
        END IF
        IF (RANK.eq.NPROCS-1) THEN
! If we have the highest plane, apply the boundary conditions
          CALL APPLY_BC_3_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C,K)
        END IF
      RETURN
      END


      SUBROUTINE APPLY_BC_REM_DIV_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INTEGER I,K
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
C Apply the boundary conditions
        IF (RANK.EQ.0) THEN
        DO I=0,NKX
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
          IF ((K.EQ.0).AND.(I.EQ.0)) THEN
C Otherwise the matrix will be singular
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
            MATL_C(I,1)=0.
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)
          ELSE
C Use Dirichlet boundary conditions, dp/dz=0 at walls
            MATL_C(I,1)=0.
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)
          END IF
        END DO
        END IF
C Apply the boundary conditions
        IF (RANK.EQ.NPROCS-1) THEN
          DO I=0,NKX
            MATL_C(I,NY)=1.
            MATD_C(I,NY)=-1.
            MATU_C(I,NY)=0.
            VEC_C(I,NY)=(0.,0.)
          END DO
        END IF


        RETURN
        END
      SUBROUTINE APPLY_BC_POISSON_MPI(MATL_C,MATD_C,MATU_C,VEC_C,K)
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'
      INTEGER I,K
! We first need to check to see which processor we are, if we are
! the upper or lowermost process, then apply boundary conditions
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
        IF (RANK.eq.0) THEN
        DO I=0,NKX
          IF ((I.EQ.0).AND.(K.EQ.0)) THEN
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)
          ELSE
! Here, apply Neumann boundary conditions (dp/dz=0) at the walls
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)
          END IF
        END DO
        END IF
C Use dirichlet boundary condition at the lower wall to
C prevent the tridiagonal matrix from becomming singular for i,k=0
        IF (RANK.eq.NPROCS-1) THEN
        DO I=0,NKX
            MATD_C(I,NY)=-1.
            MATL_C(I,NY)=1.
            VEC_C(I,NY)=(0.,0.)
        END DO
        END IF

      RETURN
      END

      SUBROUTINE APPLY_BC_VEL_MPI
! This subroutine applies the boundary conditions for the Poisson Eq.
! Note, MATL, MATD, etc. are dimensioned in header
      INCLUDE 'header'

! Apply Boundary conditions to velocity field
      IF (RANK.EQ.0) THEN
        CALL APPLY_BC_VEL_LOWER
      END IF
      IF (RANK.EQ.NPROCS-1) THEN
        CALL APPLY_BC_VEL_UPPER
      END IF

      RETURN
      END


      SUBROUTINE MPI_IO_READ(FNAME,FNAME_TH,FNAME_P)
      INCLUDE 'header'

      CHARACTER*60 FNAME
      CHARACTER*60 FNAME_P
      CHARACTER*60 FNAME_TH(N_TH)
      INTEGER I, J, K, N
      INTEGER NUM_READ_X

      integer linenum

      INTEGER fh,fh_th(N_TH),fh_p

      complex*16 test(0:100)

! Read in velocity field from diablo.saved
      IF (NUM_PER_DIR.EQ.3) THEN
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD, FNAME, MPI_MODE_RDONLY,
     &                   MPI_INFO_NULL, fh, IERROR)
! First, read in the header
      offset=4
      CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_SET,IERROR)
      CALL MPI_FILE_READ(fh, NX_T, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_READ(fh, NY_T, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_READ(fh, NZ_T, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_READ(fh, NUM_PER_DIR_T, 1, MPI_INTEGER
     &                   ,status, IERROR)
      CALL MPI_FILE_READ(fh, TIME, 1, MPI_DOUBLE_PRECISION
     &                   , status, IERROR)
      CALL MPI_FILE_READ(fh, TIME_STEP, 1, MPI_INTEGER
     &                   , status, IERROR)
! Check to make sure that the array dimensions match
      NKX_T=NX_T/3
      TNKZ_T=2*(NZ_T/3)
      IF ((NX .NE. NX_T).OR.(NY .NE. NY_T).OR.(NZ .NE. NZ_T)) then
           if (RANK.eq.0) write(*,*) 'NX,NY,NZ: ',NX,NY,NZ
           if (RANK.eq.0) write(*,*) 'NX_T,NY_T,NZ_T: ',NX_T,NY_T,NZ_T
           STOP 'Error: old flowfield wrong dimensions. '
      END IF
      IF (NUM_PER_DIR .NE. NUM_PER_DIR_T)
     *     STOP 'Error: old flowfield wrong NUM_PER_DIR. '

! The number of values to read in over x for our rank
      NUM_READ_X=MIN(NKX_S+1,(NKX+1-MOD(RANK,NP_S)*(NKX_S+1)))
! Skip the following offset
       offset=8
       CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)

       linenum=1

! Read in CU1
       do j=0,TNKY
! First set the initial offset in z for our RANK
        OFFSET=((RANK/NP_S)*(NKX+1)*(TNKZ_S+1))*16
        linenum=linenum+((RANK/NP_S)*(NKX+1)*(TNKZ_S+1))
        CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
        k=0
        do while ((k.le.TNKZ_S)
     &           .AND.((RANK/NP_S)*(TNKZ_S+1)+K.le.TNKZ))
! Skip the beginning of the row in x
          OFFSET=(MOD(RANK,NP_S)*(NKX_S+1))*16
        linenum=linenum+(MOD(RANK,NP_S)*(NKX_S+1))
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
          CALL MPI_FILE_READ(fh,CU1(0:NUM_READ_X-1,k,j),NUM_READ_X
     &                    ,MPI_DOUBLE_COMPLEX,status,IERROR)
        linenum=linenum+NUM_READ_X
! Skip the rest of the row in x
          OFFSET=((NKX+1)-NUM_READ_X-(NKX_S+1)*MOD(RANK,NP_S))*16
          linenum=linenum+((NKX+1)-NUM_READ_X-(NKX_S+1)*MOD(RANK,NP_S))
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
          k=k+1
        end do
! If we aren't at the top of the column in k
        if ((RANK/NP_S)*(TNKZ_S+1)+K-1.lt.TNKZ) then
! Skip the rest of this j index
          OFFSET=(TNKZ-((RANK/NP_S)*(TNKZ_S+1)+K-1))*(NKX+1)*16
          linenum=linenum+(TNKZ-((RANK/NP_S)*(TNKZ_S+1)+K-1))*(NKX+1)
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
        end if
      end do

! Read in CU2
       do j=0,TNKY
! First set the initial offset in z for our RANK
        OFFSET=((RANK/NP_S)*(NKX+1)*(TNKZ_S+1))*16
        CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
        k=0
        do while ((k.le.TNKZ_S)
     &           .AND.((RANK/NP_S)*(TNKZ_S+1)+K.le.TNKZ))
! Skip the beginning of the row in x
          OFFSET=(MOD(RANK,NP_S)*(NKX_S+1))*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
          CALL MPI_FILE_READ(fh,CU2(0:NUM_READ_X-1,k,j),NUM_READ_X
     &                    ,MPI_DOUBLE_COMPLEX,status,IERROR)
! Skip the rest of the row in x
          OFFSET=((NKX+1)-NUM_READ_X-(NKX_S+1)*MOD(RANK,NP_S))*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
          k=k+1
        end do
! If we aren't at the top of the column in k
        if ((RANK/NP_S)*(TNKZ_S+1)+K-1.lt.TNKZ) then
! Skip the rest of this j index
          OFFSET=(TNKZ-((RANK/NP_S)*(TNKZ_S+1)+K-1))*(NKX+1)*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
        end if
      end do

! Read in CU3
       do j=0,TNKY
! First set the initial offset in z for our RANK
        OFFSET=((RANK/NP_S)*(NKX+1)*(TNKZ_S+1))*16
        CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
        k=0
        do while ((k.le.TNKZ_S)
     &           .AND.((RANK/NP_S)*(TNKZ_S+1)+K.le.TNKZ))
! Skip the beginning of the row in x
          OFFSET=(MOD(RANK,NP_S)*(NKX_S+1))*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
          CALL MPI_FILE_READ(fh,CU3(0:NUM_READ_X-1,k,j),NUM_READ_X
     &                    ,MPI_DOUBLE_COMPLEX,status,IERROR)
! Skip the rest of the row in x
          OFFSET=((NKX+1)-NUM_READ_X-(NKX_S+1)*MOD(RANK,NP_S))*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
          k=k+1
        end do
! If we aren't at the top of the column in k
        if ((RANK/NP_S)*(TNKZ_S+1)+K-1.lt.TNKZ) then
! Skip the rest of this j index
          OFFSET=(TNKZ-((RANK/NP_S)*(TNKZ_S+1)+K-1))*(NKX+1)*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
        end if
      end do




! Now, read in the scalar files, diablo_th01.start, etc.

      DO N=1,NUM_READ_TH
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD, FNAME_TH(N), MPI_MODE_RDONLY,
     &                   MPI_INFO_NULL, fh_th(n), IERROR)
! First, read in the header
      offset=4
      CALL MPI_FILE_SEEK(fh_th(n),OFFSET,MPI_SEEK_SET,IERROR)
      CALL MPI_FILE_READ(fh_th(n), NX_T, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_READ(fh_th(n), NY_T, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_READ(fh_th(n), NZ_T, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_READ(fh_th(n), NUM_PER_DIR_T, 1, MPI_INTEGER
     &                   ,status, IERROR)
      CALL MPI_FILE_READ(fh_th(n), TIME, 1, MPI_DOUBLE_PRECISION
     &                   , status, IERROR)
      CALL MPI_FILE_READ(fh_th(n), TIME_STEP, 1, MPI_INTEGER
     &                   , status, IERROR)

! The number of values to read in over x for our rank
      NUM_READ_X=MIN(NKX_S_TH+1,(NKX_TH+1-MOD(RANK,NP_S)*(NKX_S_TH+1)))
! Skip the following offset
       offset=8
       CALL MPI_FILE_SEEK(fh_th(n),OFFSET,MPI_SEEK_CUR,IERROR)
! Read in CTH
       do j=0,TNKY_TH
! First set the initial offset in z for our RANK
        OFFSET=((RANK/NP_S)*(NKX_TH+1)*(TNKZ_S_TH+1))*16
        CALL MPI_FILE_SEEK(fh_th(n),OFFSET,MPI_SEEK_CUR,IERROR)
        k=0
        do while ((k.le.TNKZ_S_TH)
     &           .AND.((RANK/NP_S)*(TNKZ_S_TH+1)+K.le.TNKZ_TH))
! Skip the beginning of the row in x
          OFFSET=(MOD(RANK,NP_S)*(NKX_S_TH+1))*16
          CALL MPI_FILE_SEEK(fh_th(n),OFFSET,MPI_SEEK_CUR,IERROR)
          CALL MPI_FILE_READ(fh_th(n),CTH(0:NUM_READ_X-1,k,j,n),
     &              NUM_READ_X,MPI_DOUBLE_COMPLEX,status,IERROR)
! Skip the rest of the row in x
          OFFSET=((NKX_TH+1)-NUM_READ_X-(NKX_S_TH+1)*MOD(RANK,NP_S))*16
          CALL MPI_FILE_SEEK(fh_th(n),OFFSET,MPI_SEEK_CUR,IERROR)
          k=k+1
        end do
! If we aren't at the top of the column in k
        if ((RANK/NP_S)*(TNKZ_S_TH+1)+K-1.lt.TNKZ_TH) then
! Skip the rest of this j index
          OFFSET=(TNKZ_TH-((RANK/NP_S)*(TNKZ_S_TH+1)+K-1))*(NKX_TH+1)*16
          CALL MPI_FILE_SEEK(fh_th(n),OFFSET,MPI_SEEK_CUR,IERROR)
        end if
      end do

! End do number of scalars
      END DO

! Now, Optionally read in the pressure field
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD, FNAME_P, MPI_MODE_RDONLY,
     &                   MPI_INFO_NULL, fh_p, IERROR)
! First, read in the header
      offset=4
      CALL MPI_FILE_SEEK(fh_p,OFFSET,MPI_SEEK_SET,IERROR)
      CALL MPI_FILE_READ(fh_p, NX_T, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_READ(fh_p, NY_T, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_READ(fh_p, NZ_T, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_READ(fh_p, NUM_PER_DIR_T, 1, MPI_INTEGER
     &                   ,status, IERROR)
      CALL MPI_FILE_READ(fh_p, TIME, 1, MPI_DOUBLE_PRECISION
     &                   , status, IERROR)
      CALL MPI_FILE_READ(fh_p, TIME_STEP, 1, MPI_INTEGER
     &                   , status, IERROR)

! The number of values to read in over x for our rank
      NUM_READ_X=MIN(NKX_S+1,(NKX+1-MOD(RANK,NP_S)*(NKX_S+1)))
! Skip the following offset
       offset=8
       CALL MPI_FILE_SEEK(fh_p,OFFSET,MPI_SEEK_CUR,IERROR)
! Read in CP
       do j=0,TNKY
! First set the initial offset in z for our RANK
        OFFSET=((RANK/NP_S)*(NKX+1)*(TNKZ_S+1))*16
        CALL MPI_FILE_SEEK(fh_p,OFFSET,MPI_SEEK_CUR,IERROR)
        k=0
        do while ((k.le.TNKZ_S)
     &           .AND.((RANK/NP_S)*(TNKZ_S+1)+K.le.TNKZ))
! Skip the beginning of the row in x
          OFFSET=(MOD(RANK,NP_S)*(NKX_S+1))*16
          CALL MPI_FILE_SEEK(fh_p,OFFSET,MPI_SEEK_CUR,IERROR)
          CALL MPI_FILE_READ(fh_p,CP(0:NUM_READ_X-1,k,j),
     &              NUM_READ_X,MPI_DOUBLE_COMPLEX,status,IERROR)
! Skip the rest of the row in x
          OFFSET=((NKX+1)-NUM_READ_X-(NKX_S+1)*MOD(RANK,NP_S))*16
          CALL MPI_FILE_SEEK(fh_p,OFFSET,MPI_SEEK_CUR,IERROR)
          k=k+1
        end do
! If we aren't at the top of the column in k
        if ((RANK/NP_S)*(TNKZ_S+1)+K-1.lt.TNKZ) then
! Skip the rest of this j index
          OFFSET=(TNKZ-((RANK/NP_S)*(TNKZ_S+1)+K-1))*(NKX+1)*16
          CALL MPI_FILE_SEEK(fh_p,OFFSET,MPI_SEEK_CUR,IERROR)
        end if
      end do

      ELSE
        write(*,*) 'ERROR: MPI-IO Not defined for NUM_PER_DIR.NE.3'
        pause
      END IF


      RETURN
      END


      SUBROUTINE MPI_IO_WRITE(FNAME,FNAME_TH,FNAME_P)
      INCLUDE 'header'

      CHARACTER*60 FNAME
      CHARACTER*60 FNAME_P
      CHARACTER*60 FNAME_TH(N_TH)
      INTEGER I, J, K, N
      INTEGER NUM_WRITE_X

      INTEGER fh,fh_th(N_TH),fh_p

      INTEGER HEADER

! Write velocity field to diablo.saved
      IF (NUM_PER_DIR.EQ.3) THEN
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD, FNAME
     &                   , MPI_MODE_CREATE + MPI_MODE_EXCL
     &                     + MPI_MODE_WRONLY
     &                   , MPI_INFO_NULL, fh, IERROR)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      IF (IERROR.NE.0) THEN
        WRITE(*,*) '****ERROR***** OUTPUT FILE EXISTS: ',FNAME
        STOP
!        write(*,*) 'WARNING, FILE EXISTS, RANK=',RANK
!        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!        IF (RANK.eq.0) THEN
!          CALL SYSTEM('rm -f '//FNAME)
!        END IF
!        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!        CALL MPI_FILE_OPEN(MPI_COMM_WORLD, FNAME
!     &                   , MPI_MODE_CREATE + MPI_MODE_EXCL
!     &                     + MPI_MODE_WRONLY
!     &                   , MPI_INFO_NULL, fh, IERROR)
      END IF

!      IF (IERROR.NE.0) THEN
!        write(*,*) 'FNAME: ',FNAME
!        PAUSE 'WARNING: IERROR.NE.0 ON MPI_FILE_OPEN'
!      END IF

! First, write the Fortran pre-header containing the number of bytes in this line
      IF (RANK.eq.0) THEN
        HEADER=28
        CALL MPI_FILE_WRITE(fh,HEADER,1,MPI_INTEGER,status,ierror)
      ELSE
! ... or skip the Fortran pre-header
        offset=4
        CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_SET,IERROR)
      END IF
      CALL MPI_FILE_WRITE(fh, NX, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_WRITE(fh, NY, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_WRITE(fh, NZ, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_WRITE(fh, NUM_PER_DIR, 1, MPI_INTEGER
     &                   ,status, IERROR)
      CALL MPI_FILE_WRITE(fh, TIME, 1, MPI_DOUBLE_PRECISION
     &                   , status, IERROR)
      CALL MPI_FILE_WRITE(fh, TIME_STEP, 1, MPI_INTEGER
     &                   , status, IERROR)
! Write the Fortran post-header containing the number of bytes in this line
      IF (RANK.eq.0) THEN
        HEADER=28
        CALL MPI_FILE_WRITE(fh,HEADER,1,MPI_INTEGER,status,ierror)
      ELSE
! ...or skip the fortran post-header
        offset=4
        CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
      END IF

! Write the Fortran pre-header containing the number of bytes
       NUM_WRITE_X=MIN(NKX_S+1,(NKX+1-MOD(RANK,NP_S)*(NKX_S+1)))

! Write the Fortran pre-header
       IF (RANK.eq.0) THEN
         HEADER=(TNKY+1)*(TNKZ+1)*(NKX+1)*16*3
         CALL MPI_FILE_WRITE(fh,HEADER,1,MPI_INTEGER,status,IERROR)
       ELSE
! ...or skip the fortran pre-header
         offset=4
         CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
       END IF

! Write CU1
       do j=0,TNKY
! First set the initial offset in z for our RANK
        OFFSET=((RANK/NP_S)*(NKX+1)*(TNKZ_S+1))*16
        CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
        k=0
        do while ((k.le.TNKZ_S)
     &           .AND.((RANK/NP_S)*(TNKZ_S+1)+K.le.TNKZ))
! Skip the beginning of the row in x
          OFFSET=(MOD(RANK,NP_S)*(NKX_S+1))*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
          CALL MPI_FILE_WRITE(fh,CU1(0:NUM_WRITE_X-1,k,j),NUM_WRITE_X
     &                    ,MPI_DOUBLE_COMPLEX,status,IERROR)
! Skip the rest of the row in x
          OFFSET=((NKX+1)-NUM_WRITE_X-(NKX_S+1)*MOD(RANK,NP_S))*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
          k=k+1
        end do
! If we aren't at the top of the column in k
        if ((RANK/NP_S)*(TNKZ_S+1)+K-1.lt.TNKZ) then
! Skip the rest of this j index
          OFFSET=(TNKZ-((RANK/NP_S)*(TNKZ_S+1)+K-1))*(NKX+1)*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
        end if
      end do

! Write CU2
       do j=0,TNKY
! First set the initial offset in z for our RANK
        OFFSET=((RANK/NP_S)*(NKX+1)*(TNKZ_S+1))*16
        CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
        k=0
        do while ((k.le.TNKZ_S)
     &           .AND.((RANK/NP_S)*(TNKZ_S+1)+K.le.TNKZ))
! Skip the beginning of the row in x
          OFFSET=(MOD(RANK,NP_S)*(NKX_S+1))*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
          CALL MPI_FILE_WRITE(fh,CU2(0:NUM_WRITE_X-1,k,j),NUM_WRITE_X
     &                    ,MPI_DOUBLE_COMPLEX,status,IERROR)
! Skip the rest of the row in x
          OFFSET=((NKX+1)-NUM_WRITE_X-(NKX_S+1)*MOD(RANK,NP_S))*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
          k=k+1
        end do
! If we aren't at the top of the column in k
        if ((RANK/NP_S)*(TNKZ_S+1)+K-1.lt.TNKZ) then
! Skip the rest of this j index
          OFFSET=(TNKZ-((RANK/NP_S)*(TNKZ_S+1)+K-1))*(NKX+1)*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
        end if
      end do

! Write CU3
       do j=0,TNKY
! First set the initial offset in z for our RANK
        OFFSET=((RANK/NP_S)*(NKX+1)*(TNKZ_S+1))*16
        CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
        k=0
        do while ((k.le.TNKZ_S)
     &           .AND.((RANK/NP_S)*(TNKZ_S+1)+K.le.TNKZ))
! Skip the beginning of the row in x
          OFFSET=(MOD(RANK,NP_S)*(NKX_S+1))*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
          CALL MPI_FILE_WRITE(fh,CU3(0:NUM_WRITE_X-1,k,j),NUM_WRITE_X
     &                    ,MPI_DOUBLE_COMPLEX,status,IERROR)
! Skip the rest of the row in x
          OFFSET=((NKX+1)-NUM_WRITE_X-(NKX_S+1)*MOD(RANK,NP_S))*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
          k=k+1
        end do
! If we aren't at the top of the column in k
        if ((RANK/NP_S)*(TNKZ_S+1)+K-1.lt.TNKZ) then
! Skip the rest of this j index
          OFFSET=(TNKZ-((RANK/NP_S)*(TNKZ_S+1)+K-1))*(NKX+1)*16
          CALL MPI_FILE_SEEK(fh,OFFSET,MPI_SEEK_CUR,IERROR)
        end if
      end do

! Make sure that all processes are done reading
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
! Write the fortran post-header
       IF (RANK.eq.0) THEN
         HEADER=(TNKY+1)*(TNKZ+1)*(NKX+1)*16*3
         CALL MPI_FILE_WRITE(fh,HEADER,1,MPI_INTEGER,status,IERROR)
       END IF


! Write to the scalar files, diablo_th01.saved, etc.

      DO N=1,N_TH
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD, FNAME_TH(N)
     &                   , MPI_MODE_CREATE + MPI_MODE_EXCL
     &                     + MPI_MODE_WRONLY
     &                   , MPI_INFO_NULL, fh_th(N), IERROR)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      IF (IERROR.NE.0) THEN
        WRITE(*,*) '****ERROR***** OUTPUT FILE EXISTS: ',FNAME_TH(N)
        STOP
!        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!        IF (RANK.eq.0) THEN
!          CALL SYSTEM('rm -f '//FNAME_TH(N))
!        END IF
!        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!        CALL MPI_FILE_OPEN(MPI_COMM_WORLD, FNAME_TH(N)
!     &                   , MPI_MODE_CREATE + MPI_MODE_EXCL
!     &                     + MPI_MODE_WRONLY
!     &                   , MPI_INFO_NULL, fh_th(N), IERROR)
      END IF

!      IF (IERROR.NE.0) THEN
!        write(*,*) 'FNAME_TH: ',FNAME_TH(N)
!        PAUSE 'WARNING: IERROR.NE.0 ON MPI_FILE_OPEN'
!      END IF

! Write the Fortran post-header containing the number of bytes in this line
      IF (RANK.eq.0) THEN
        HEADER=28
        CALL MPI_FILE_WRITE(fh_th(n),HEADER,1,MPI_INTEGER,status,ierror)
      ELSE
! ...or skip the Fortran post-header
        offset=4
        CALL MPI_FILE_SEEK(fh_th(n),OFFSET,MPI_SEEK_SET,IERROR)
      END IF

      CALL MPI_FILE_WRITE(fh_th(N),NX_TH,1,MPI_INTEGER,status,IERROR)
      CALL MPI_FILE_WRITE(fh_th(N),NY_TH,1,MPI_INTEGER,status,IERROR)
      CALL MPI_FILE_WRITE(fh_th(N),NZ_TH,1,MPI_INTEGER,status,IERROR)
      CALL MPI_FILE_WRITE(fh_th(N),NUM_PER_DIR, 1, MPI_INTEGER
     &                   ,status, IERROR)
      CALL MPI_FILE_WRITE(fh_th(N), TIME, 1, MPI_DOUBLE_PRECISION
     &                   , status, IERROR)
      CALL MPI_FILE_WRITE(fh_th(N), TIME_STEP, 1, MPI_INTEGER
     &                   , status, IERROR)
! Write the Fortran post-header containing the number of bytes in this line
      IF (RANK.eq.0) THEN
        HEADER=28
        CALL MPI_FILE_WRITE(fh_th(n),HEADER,1,MPI_INTEGER,status,ierror)
      ELSE
! ...or skip the Fortran post-header
        offset=4
        CALL MPI_FILE_SEEK(fh_th(n),OFFSET,MPI_SEEK_CUR,IERROR)
      END IF

! The number of values to write in x for our rank
      NUM_WRITE_X=MIN(NKX_S_TH+1,(NKX_TH+1-MOD(RANK,NP_S)*(NKX_S_TH+1)))
! Write the Fortran pre-header
       IF (RANK.eq.0) THEN
         HEADER=(TNKY_TH+1)*(TNKZ_TH+1)*(NKX_TH+1)*16
         CALL MPI_FILE_WRITE(fh_th(n),HEADER,1
     &                       ,MPI_INTEGER,status,IERROR)
       ELSE
! ...or skip the Fortran pre-header
         offset=4
         CALL MPI_FILE_SEEK(fh_th(n),OFFSET,MPI_SEEK_CUR,IERROR)
       END IF

! Write CTH
       do j=0,TNKY_TH
! First set the initial offset in z for our RANK
        OFFSET=((RANK/NP_S)*(NKX_TH+1)*(TNKZ_S_TH+1))*16
        CALL MPI_FILE_SEEK(fh_th(n),OFFSET,MPI_SEEK_CUR,IERROR)
        k=0
        do while ((k.le.TNKZ_S_TH)
     &           .AND.((RANK/NP_S)*(TNKZ_S_TH+1)+K.le.TNKZ_TH))
! Skip the beginning of the row in x
          OFFSET=(MOD(RANK,NP_S)*(NKX_S_TH+1))*16
          CALL MPI_FILE_SEEK(fh_th(n),OFFSET,MPI_SEEK_CUR,IERROR)
          CALL MPI_FILE_WRITE(fh_th(n),CTH(0:NUM_WRITE_X-1,k,j,n),
     &              NUM_WRITE_X,MPI_DOUBLE_COMPLEX,status,IERROR)
! Skip the rest of the row in x
          OFFSET=((NKX_TH+1)-NUM_WRITE_X-(NKX_S_TH+1)*MOD(RANK,NP_S))*16
          CALL MPI_FILE_SEEK(fh_th(n),OFFSET,MPI_SEEK_CUR,IERROR)
          k=k+1
        end do
! If we aren't at the top of the column in k
        if ((RANK/NP_S)*(TNKZ_S_TH+1)+K-1.lt.TNKZ_TH) then
! Skip the rest of this j index
          OFFSET=(TNKZ_TH-((RANK/NP_S)*(TNKZ_S_TH+1)+K-1))*(NKX_TH+1)*16
          CALL MPI_FILE_SEEK(fh_th(n),OFFSET,MPI_SEEK_CUR,IERROR)
        end if
      end do
! Make sure that all processes are done reading
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
! Write the Fortran post-header
       IF (RANK.eq.0) THEN
         HEADER=(TNKY_TH+1)*(TNKZ_TH+1)*(NKX_TH+1)*16
         CALL MPI_FILE_WRITE(fh_th(n),HEADER,1
     &                       ,MPI_INTEGER,status,IERROR)
       END IF
! End do number of scalars
      END DO

! Optionally, write out the pressure field

      CALL MPI_FILE_OPEN(MPI_COMM_WORLD, FNAME_P
     &                   , MPI_MODE_CREATE + MPI_MODE_EXCL
     &                     + MPI_MODE_WRONLY
     &                   , MPI_INFO_NULL, fh_p, IERROR)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      IF (IERROR.NE.0) THEN
        WRITE(*,*) '****ERROR***** OUTPUT P FILE EXISTS: ',FNAME_P
        STOP
!        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!        IF (RANK.eq.0) THEN
!          CALL SYSTEM('rm -f '//FNAME_P)
!        END IF
!        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!        CALL MPI_FILE_OPEN(MPI_COMM_WORLD, FNAME_P
!     &                   , MPI_MODE_CREATE + MPI_MODE_EXCL
!     &                     + MPI_MODE_WRONLY
!     &                   , MPI_INFO_NULL, fh_p, IERROR)
      END IF
!      IF (IERROR.NE.0) THEN
!        write(*,*) 'FNAME_P: ',FNAME_P
!        PAUSE 'WARNING: IERROR.NE.0 ON MPI_FILE_OPEN'
!      END IF
! Write the Fortran post-header containing the number of bytes in this line
      IF (RANK.eq.0) THEN
        HEADER=28
        CALL MPI_FILE_WRITE(fh_p,HEADER,1,MPI_INTEGER,status,ierror)
      ELSE
! ...or skip the Fortran post-header
        offset=4
        CALL MPI_FILE_SEEK(fh_p,OFFSET,MPI_SEEK_SET,IERROR)
      END IF

      CALL MPI_FILE_WRITE(fh_p, NX, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_WRITE(fh_p, NY, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_WRITE(fh_p, NZ, 1, MPI_INTEGER, status, IERROR)
      CALL MPI_FILE_WRITE(fh_p, NUM_PER_DIR, 1, MPI_INTEGER
     &                   ,status, IERROR)
      CALL MPI_FILE_WRITE(fh_p, TIME, 1, MPI_DOUBLE_PRECISION
     &                   , status, IERROR)
      CALL MPI_FILE_WRITE(fh_p, TIME_STEP, 1, MPI_INTEGER
     &                   , status, IERROR)
! Write the Fortran post-header containing the number of bytes in this line
      IF (RANK.eq.0) THEN
        HEADER=28
        CALL MPI_FILE_WRITE(fh_p,HEADER,1,MPI_INTEGER,status,ierror)
      ELSE
! ...or skip the Fortran post-header
        offset=4
        CALL MPI_FILE_SEEK(fh_p,OFFSET,MPI_SEEK_CUR,IERROR)
      END IF
! The number of values to write in x for our rank
      NUM_WRITE_X=MIN(NKX_S+1,(NKX+1-MOD(RANK,NP_S)*(NKX_S+1)))
! Write the Fortran pre-header
       IF (RANK.eq.0) THEN
         HEADER=(TNKY+1)*(TNKZ+1)*(NKX+1)*16
         CALL MPI_FILE_WRITE(fh_p,HEADER,1
     &                       ,MPI_INTEGER,status,IERROR)
       ELSE
! ...or skip the Fortran pre-header
         offset=4
         CALL MPI_FILE_SEEK(fh_p,OFFSET,MPI_SEEK_CUR,IERROR)
       END IF

! Write CP
       do j=0,TNKY
! First set the initial offset in z for our RANK
        OFFSET=((RANK/NP_S)*(NKX+1)*(TNKZ_S+1))*16
        CALL MPI_FILE_SEEK(fh_p,OFFSET,MPI_SEEK_CUR,IERROR)
        k=0
        do while ((k.le.TNKZ_S)
     &           .AND.((RANK/NP_S)*(TNKZ_S+1)+K.le.TNKZ))
! Skip the beginning of the row in x
          OFFSET=(MOD(RANK,NP_S)*(NKX_S+1))*16
          CALL MPI_FILE_SEEK(fh_p,OFFSET,MPI_SEEK_CUR,IERROR)
          CALL MPI_FILE_WRITE(fh_p,CP(0:NUM_WRITE_X-1,k,j),
     &              NUM_WRITE_X,MPI_DOUBLE_COMPLEX,status,IERROR)
! Skip the rest of the row in x
          OFFSET=((NKX+1)-NUM_WRITE_X-(NKX_S+1)*MOD(RANK,NP_S))*16
          CALL MPI_FILE_SEEK(fh_p,OFFSET,MPI_SEEK_CUR,IERROR)
          k=k+1
        end do
! If we aren't at the top of the column in k
        if ((RANK/NP_S)*(TNKZ_S+1)+K-1.lt.TNKZ) then
! Skip the rest of this j index
          OFFSET=(TNKZ-((RANK/NP_S)*(TNKZ_S+1)+K-1))*(NKX+1)*16
          CALL MPI_FILE_SEEK(fh_p,OFFSET,MPI_SEEK_CUR,IERROR)
        end if
      end do
! Make sure that all processes are done reading
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
! Write the Fortran post-header
       IF (RANK.eq.0) THEN
         HEADER=(TNKY+1)*(TNKZ+1)*(NKX+1)*16
         CALL MPI_FILE_WRITE(fh_p,HEADER,1
     &                       ,MPI_INTEGER,status,IERROR)
       END IF

      ELSE
        write(*,*) 'ERROR: MPI-IO Not defined for NUM_PER_DIR.NE.3'
        pause
      END IF

      RETURN
      END
