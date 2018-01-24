

!----*|--.---------.---------.---------.---------.---------.---------.-|------
      SUBROUTINE INIT_MPI_TH
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      INCLUDE 'header'
      INCLUDE 'mpif.h'
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables

        NXM_S_TH=NX_S_TH-1
        NYM_S_TH=NY_S_TH-1
        NZM_S_TH=NZ_S_TH-1

      RETURN
      END


      SUBROUTINE TRANSPOSE_MPI_Z_TO_X_TH(CA_Z,CA_X)
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables


      complex*16 CA_Z(0:NX_S_TH/2,0:NZ_TH+1,0:NY_S_TH)
      complex*16 CA_X(0:NX_TH/2,0:NZ_S_TH,0:NY_S_TH)

      complex*16 SENDBUF(0:NKX_S_TH,0:NZ_S_TH,0:NY_S_TH)
      complex*16 RECVBUF(0:NKX_S_TH,0:NZ_S_TH,0:NY_S_TH)

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

      sendcount=(NKX_S_TH+1)*(NZ_S_TH+1)*(NY_S_TH+1)
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
        K1=(NZ_S_TH+1)*BLOCK
        K2=NZ_TH+1-K1
        DO J=0,NY_S_TH
          DO K=0,MIN(K2,NZ_S_TH)
            DO I=0,NKX_S_TH
              SENDBUF(I,K,J)=CA_Z(I,K+K1,J)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself
          RECVBUF=0.d0
          CALL MPI_SENDRECV(
     &        sendbuf(0,0,0),sendcount,MPI_DOUBLE_COMPLEX
     &               ,block+(RANK/NP_S)*NP_S,12
     &        ,recvbuf(0,0,0),recvcount,MPI_DOUBLE_COMPLEX
     &               ,block+(RANK/NP_S)*NP_S,12
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        I1=(NKX_S_TH+1)*BLOCK
        I2=NKX_S_TH+I1
        DO J=0,NY_S_TH
          DO K=0,NZ_S_TH
            DO I=I1,MIN(I2,NX_TH/2)
! Make sure that we don't write beyond the end of the array
              CA_X(I,K,J)=RECVBUF(I-I1,K,J)
            END DO
          END DO
        END DO

        ELSE
! Here block=diag, we just need to copy locally from A->B
          IF (RANK.eq.0) THEN
! If we are on RANK=0, just copy directly from input to output
            DO J=0,NY_S_TH
              DO K=0,NZ_S_TH
                DO I=0,NKX_S_TH
                  CA_X(I,K,J)=CA_Z(I,K,J)
                END DO
              END DO
            END DO
          ELSE
! Here, RANK.ne.0, so we need to shift the input and output
            K1=(NZ_S_TH+1)*BLOCK
            K2=NZ_TH+1-K1
            I1=(NKX_S_TH+1)*BLOCK
            I2=I1+NKX_S_TH
            DO J=0,NY_S_TH
              DO K=0,MIN(K2,NZ_S_TH)
                DO I=I1,MIN(I2,NX_TH/2)
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

      SUBROUTINE TRANSPOSE_MPI_X_TO_Z_TH(CA_X,CA_Z)
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables

      COMPLEX*16 CA_Z(0:NKX_S_TH,0:NZ_TH+1,0:NY_S_TH)
      COMPLEX*16 CA_X(0:NX_TH/2,0:NZ_S_TH,0:NY_S_TH)

      COMPLEX*16 SENDBUF(0:NKX_S_TH,0:NZ_S_TH,0:NY_S_TH)
      COMPLEX*16 RECVBUF(0:NKX_S_TH,0:NZ_S_TH,0:NY_S_TH)

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

      sendcount=(NKX_S_TH+1)*(NZ_S_TH+1)*(NY_S_TH+1)
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
        I1=(NKX_S_TH+1)*BLOCK
        I2=NX_TH/2-I1
        DO J=0,NY_S_TH
          DO K=0,NZ_S_TH
            DO I=0,MIN(I2,NKX_S_TH)
              SENDBUF(I,K,J)=CA_X(I+I1,K,J)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself
          RECVBUF=0.d0
          CALL MPI_SENDRECV(
     &        sendbuf(0,0,0),sendcount,MPI_DOUBLE_COMPLEX
     &               ,block+(RANK/NP_S)*NP_S,13
     &        ,recvbuf(0,0,0),recvcount,MPI_DOUBLE_COMPLEX
     &               ,block+(RANK/NP_S)*NP_S,13
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        K1=(NZ_S_TH+1)*BLOCK
        K2=NZ_S_TH+K1
        DO J=0,NY_S_TH
          DO K=K1,MIN(K2,NZ_TH+1)
            DO I=0,NKX_S_TH
! Make sure that we don't write beyond the end of the array
              CA_Z(I,K,J)=RECVBUF(I,K-K1,J)
            END DO
          END DO
        END DO

      ELSE
! Here block=diag, we just need to copy locally from A->B
        IF (RANK.eq.0) THEN
          DO J=0,NY_S_TH
            DO K=0,NZ_S_TH
              DO I=0,NKX_S_TH
                CA_Z(I,K,J)=CA_X(I,K,J)
              END DO
            END DO
          END DO
        ELSE
          I1=(NKX_S_TH+1)*BLOCK
          K1=(NZ_S_TH+1)*BLOCK
          K2=NZ_S_TH+K1
          DO J=0,NY_S_TH
            DO K=K1,MIN(K2,NZ_TH+1)
              DO I=0,NKX_S_TH
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

      SUBROUTINE TRANSPOSE_MPI_Y_TO_Z_TH(CA_Y,CA_Z)
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables


      complex*16 CA_Z(0:NX_S_TH/2,0:NZ_TH+1,0:NY_S_TH)
      complex*16 CA_Y(0:NX_S_TH/2,0:NZ_S_TH,0:NY_TH+1)

      complex*16 SENDBUF(0:NKX_S_TH,0:TNKZ_S_TH,0:NY_S_TH)
      complex*16 RECVBUF(0:NKX_S_TH,0:TNKZ_S_TH,0:NY_S_TH)

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

      sendcount=(NKX_S_TH+1)*(TNKZ_S_TH+1)*(NY_S_TH+1)
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
        J1=(NY_S_TH+1)*BLOCK
        J2=NY_TH+1-J1
        DO J=0,MIN(J2,NY_S_TH)
          DO K=0,TNKZ_S_TH
            DO I=0,NKX_S_TH
              SENDBUF(I,K,J)=CA_Y(I,K,J+J1)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself
          RECVBUF=(0.d0,0.d0)
          CALL MPI_SENDRECV(
     &        sendbuf,sendcount,MPI_DOUBLE_COMPLEX
     &               ,MOD(RANK,NP_S)+block*NP_S,15
     &        ,recvbuf,recvcount,MPI_DOUBLE_COMPLEX
     &               ,MOD(RANK,NP_S)+block*NP_S,15
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        K1=(TNKZ_S_TH+1)*BLOCK
        K2=TNKZ_S_TH+K1
        DO J=0,NY_S_TH
          DO K=K1,MIN(K2,NZ_TH+1)
            DO I=0,NKX_S_TH
! Make sure that we don't write beyond the end of the array
              CA_Z(I,K,J)=RECVBUF(I,K-K1,J)
            END DO
          END DO
        END DO

        ELSE
! Here block=diag, we just need to copy locally from A->B
          IF (RANK.eq.0) THEN
          DO J=0,NY_S_TH
            DO K=0,TNKZ_S_TH
              DO I=0,NKX_S_TH
                CA_Z(I,K,J)=CA_Y(I,K,J)
              END DO
            END DO
          END DO
          ELSE
          J1=(NY_S_TH+1)*BLOCK
          J2=NY_TH+1-J1
          K1=(TNKZ_S_TH+1)*BLOCK
          K2=K1+TNKZ_S_TH
          DO J=0,MIN(J2,NY_S_TH)
            DO K=K1,MIN(K2,NZ_TH+1)
              DO I=0,NKX_S_TH
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

      SUBROUTINE TRANSPOSE_MPI_Z_TO_Y_TH(CA_Z,CA_Y)
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables

      complex*16 CA_Z(0:NKX_S_TH,0:NZ_TH+1,0:NY_S_TH)
      complex*16 CA_Y(0:NKX_S_TH,0:TNKZ_S_TH,0:NY_TH+1)

      complex*16 SENDBUF(0:NKX_S_TH,0:TNKZ_S_TH,0:NY_S_TH)
      complex*16 RECVBUF(0:NKX_S_TH,0:TNKZ_S_TH,0:NY_S_TH)

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

      sendcount=(NKX_S_TH+1)*(TNKZ_S_TH+1)*(NY_S_TH+1)
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
        K1=(TNKZ_S_TH+1)*BLOCK
        K2=NZ_TH+1-K1
        DO J=0,NY_S_TH
          DO K=0,MIN(K2,TNKZ_S_TH)
            DO I=0,NKX_S_TH
              SENDBUF(I,K,J)=CA_Z(I,K+K1,J)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself
          RECVBUF=0.d0
          CALL MPI_SENDRECV(
     &        sendbuf,sendcount,MPI_DOUBLE_COMPLEX
     &               ,MOD(RANK,NP_S)+block*NP_S,14
     &        ,recvbuf,recvcount,MPI_DOUBLE_COMPLEX
     &               ,MOD(RANK,NP_S)+block*NP_S,14
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        J1=(NY_S_TH+1)*BLOCK
        J2=NY_S_TH+J1
        DO J=J1,MIN(J2,NY_TH+1)
          DO K=0,TNKZ_S_TH
            DO I=0,NKX_S_TH
! Make sure that we don't write beyond the end of the array
               CA_Y(I,K,J)=RECVBUF(I,K,J-J1)
            END DO
          END DO
        END DO

        ELSE
! Here block=diag, we just need to copy locally from A->B
          IF (RANK.eq.0) THEN
          DO J=0,NY_S_TH
              DO K=0,TNKZ_S_TH
                DO I=0,NKX_S_TH
                  CA_Y(I,K,J)=CA_Z(I,K,J)
                END DO
              END DO
          END DO
          ELSE
          K1=(TNKZ_S_TH+1)*BLOCK
          K2=NZ_TH+1-K1
          J1=(NY_S_TH+1)*BLOCK
          J2=NY_S_TH+J1
          DO J=J1,MIN(J2,NY_TH+1)
              DO K=0,MIN(K2,TNKZ_S_TH)
                DO I=0,NKX_S_TH
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



! The following subroutine is for taking a velocity field in Fourier space
! and interpolating it to the scalar (TH) grid in Physical space

      SUBROUTINE TRANSPOSE_MPI_Y_TO_Z_INTERP(CA_Y,CA_Z)
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables


      complex*16 CA_Z(0:NX_S/2,0:NZ+1,0:NY_S_TH)
      complex*16 CA_Y(0:NX_S/2,0:NZ_S,0:NY_TH+1)

      complex*16 SENDBUF(0:NKX_S,0:TNKZ_S,0:NY_S_TH)
      complex*16 RECVBUF(0:NKX_S,0:TNKZ_S,0:NY_S_TH)

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

      sendcount=(NKX_S+1)*(TNKZ_S+1)*(NY_S_TH+1)
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
        J1=(NY_S_TH+1)*BLOCK
        J2=NY_TH+1-J1
        DO J=0,MIN(J2,NY_S_TH)
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
     &               ,MOD(RANK,NP_S)+block*NP_S,15
     &        ,recvbuf,recvcount,MPI_DOUBLE_COMPLEX
     &               ,MOD(RANK,NP_S)+block*NP_S,15
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        K1=(TNKZ_S+1)*BLOCK
        K2=TNKZ_S+K1
        DO J=0,NY_S_TH
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
          DO J=0,NY_S_TH
            DO K=0,TNKZ_S
              DO I=0,NKX_S
                CA_Z(I,K,J)=CA_Y(I,K,J)
              END DO
            END DO
          END DO
          ELSE
          J1=(NY_S_TH+1)*BLOCK
          J2=NY_TH+1-J1
          K1=(TNKZ_S+1)*BLOCK
          K2=K1+TNKZ_S
          DO J=0,MIN(J2,NY_S_TH)
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


      SUBROUTINE TRANSPOSE_MPI_Z_TO_X_INTERP(CA_Z,CA_X)
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables


      complex*16 CA_Z(0:NX_S/2,0:NZ_TH+1,0:NY_S_TH)
      complex*16 CA_X(0:NX/2,0:NZ_S_TH,0:NY_S_TH)

      complex*16 SENDBUF(0:NKX_S,0:NZ_S_TH,0:NY_S_TH)
      complex*16 RECVBUF(0:NKX_S,0:NZ_S_TH,0:NY_S_TH)

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

      sendcount=(NKX_S+1)*(NZ_S_TH+1)*(NY_S_TH+1)
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
        K1=(NZ_S_TH+1)*BLOCK
        K2=NZ_TH+1-K1
        DO J=0,NY_S_TH
          DO K=0,MIN(K2,NZ_S_TH)
            DO I=0,NKX_S
              SENDBUF(I,K,J)=CA_Z(I,K+K1,J)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself
          RECVBUF=0.d0
          CALL MPI_SENDRECV(
     &        sendbuf(0,0,0),sendcount,MPI_DOUBLE_COMPLEX
     &               ,block+(RANK/NP_S)*NP_S,12
     &        ,recvbuf(0,0,0),recvcount,MPI_DOUBLE_COMPLEX
     &               ,block+(RANK/NP_S)*NP_S,12
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        I1=(NKX_S+1)*BLOCK
        I2=NKX_S+I1
        DO J=0,NY_S_TH
          DO K=0,NZ_S_TH
            DO I=I1,MIN(I2,NX/2)
! Make sure that we don't write beyond the end of the array
              CA_X(I,K,J)=RECVBUF(I-I1,K,J)
            END DO
          END DO
        END DO

        ELSE
! Here block=diag, we just need to copy locally from A->B
          IF (RANK.eq.0) THEN
! If we are on RANK=0, just copy directly from input to output
            DO J=0,NY_S_TH
              DO K=0,NZ_S_TH
                DO I=0,NKX_S
                  CA_X(I,K,J)=CA_Z(I,K,J)
                END DO
              END DO
            END DO
          ELSE
! Here, RANK.ne.0, so we need to shift the input and output
            K1=(NZ_S_TH+1)*BLOCK
            K2=NZ_TH+1-K1
            I1=(NKX_S+1)*BLOCK
            I2=I1+NKX_S
            DO J=0,NY_S_TH
              DO K=0,MIN(K2,NZ_S_TH)
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

