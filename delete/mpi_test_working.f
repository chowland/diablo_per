

      program mpi_test

      INCLUDE "/usr/local/mpich2/include/mpif.h"
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables

      integer NX,NY,NZ
      parameter (NX=16,NY=10,NZ=12)
      parameter (NX_S=3,NY_S=1,NZ_S=2)

      complex*16 A_Z(0:NX_S,0:NZ,0:NY_S)
      complex*16 A_X(0:NX,0:NZ_S,0:NY_S)
      complex*16 A_Y(0:NX_S,0:NZ_S,0:NY)

      EQUIVALENCE(A_Z,A_X,A_Y)

      CALL MPI_INIT(IERROR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERROR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,RANK,IERROR)


      do i=0,NX
        do k=0,NZ_S
          do j=0,NY_S
            A_X(i,k,j)=RANK*(NY_S+1)*(NX+1)+j+(NY_S+1)*(i)
          end do
        end do
      end do

        k=1
        do i=0,NX
        do j=0,NY_S
          write(100+rank,*) 'i,j,A_X(i,k): ',i,j,A_X(i,k,j)
        end do
        end do

      CALL TRANSPOSE_MPI_X_TO_Y(A_X,A_Y,NX,NY,NZ,NX_S,NY_S,NZ_S)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      CALL TRANSPOSE_MPI_Y_TO_X(A_Y,A_X,NX,NY,NZ,NX_S,NY_S,NZ_S)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

      k=1
        do i=0,NX
          do j=0,NY_S
            write(200+rank,*) 'i,j,A_X(i,k): ',i,j,A_X(i,k,j)
          end do
        end do


      CALL MPI_Finalize(IERROR)

      STOP
      END


      SUBROUTINE TRANSPOSE_MPI_Z_TO_X(A_Z,A_X,NX,NY,NZ,NX_S,NY_S,NZ_S)
      INCLUDE "/usr/local/mpich2/include/mpif.h"
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables

      integer NX,NY,NZ,NX_S,NY_S,NZ_S

      real*8 A_Z(0:NX_S,0:NZ,0:NY_S)
      real*8 A_X(0:NX,0:NZ_S,0:NY_S)
      real*8 B(0:NX,0:NZ_S,0:NY_S)

      real*8 SENDBUF(0:NX_S,0:NZ_S,0:NY_S)
      real*8 RECVBUF(0:NX_S,0:NZ_S,0:NY_S)

      integer DIAG
      integer block
      integer sendcount,recvcount

      NP_S=NPROCS**0.5

      DIAG=MOD(RANK,NP_S)

      sendcount=(NX_S+1)*(NZ_S+1)*(NY_S+1)
      recvcount=sendcount

      DO block=0,NP_S-1

! Get a block from the source array to be sent to another process
        DO I=0,NX_S
          DO K=0,NZ_S
            DO J=0,NY_S
              SENDBUF(I,K,J)=A_Z(I,K+(NZ_S+1)*BLOCK,J)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself


          CALL MPI_SENDRECV(
     &        sendbuf,sendcount,MPI_DOUBLE_PRECISION
     &               ,RANK+block-diag,RANK+block
     &        ,recvbuf,recvcount,MPI_DOUBLE_PRECISION
     &               ,RANK+block-diag,RANK+block
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        DO I=0,NX_S
          DO K=0,NZ_S
            DO J=0,NY_S
              IF (I+(NX_S+1)*BLOCK.le.NX) THEN
! Make sure that we don't write beyond the end of the array
                B(I+(NX_S+1)*BLOCK,K,J)=RECVBUF(I,K,J)
              END IF
            END DO
          END DO
        END DO

        ELSE
! Here block=diag, we just need to copy locally from A->B

          DO I=0,NX_S
            DO K=0,NZ_S
              DO J=0,NY_S
                sendbuf(I,K,J)=A_Z(I,K+(NZ_S+1)*BLOCK,J)
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

      SUBROUTINE TRANSPOSE_MPI_X_TO_Z(A_X,A_Z,NX,NY,NZ,NX_S,NY_S,NZ_S)
      INCLUDE "/usr/local/mpich2/include/mpif.h"
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables

      integer NX,NY,NZ,NX_S,NY_S,NZ_S

      real*8 A_Z(0:NX_S,0:NZ,0:NY_S)
      real*8 A_X(0:NX,0:NZ_S,0:NY_S)
      real*8 B(0:NX_S,0:NZ,0:NY_S)

      real*8 SENDBUF(0:NX_S,0:NZ_S,0:NY_S)
      real*8 RECVBUF(0:NX_S,0:NZ_S,0:NY_S)

      integer DIAG
      integer block
      integer sendcount,recvcount

      write(*,*) 'In X->Z'

      NP_S=NPROCS**0.5

      DIAG=MOD(RANK,NP_S)

      sendcount=(NX_S+1)*(NZ_S+1)*(NY_S+1)
      recvcount=sendcount

      DO block=0,NP_S-1


! Get a block from the source array to be sent to another process
        DO I=0,NX_S
          DO K=0,NZ_S
            DO J=0,NY_S
              SENDBUF(I,K,J)=A_X(I+(NX_S+1)*BLOCK,K,J)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself

        if ((rank.eq.0).and.(block.eq.1)) then 
          do i=0,NX_S
          do k=0,NZ_S
            write(*,*) 'i,k,sendbuf: ',i,k,sendbuf(i,k,1)
          end do
          end do
        end if

          CALL MPI_SENDRECV(
     &        sendbuf,sendcount,MPI_DOUBLE_PRECISION
     &               ,RANK+block-diag,RANK+block
     &        ,recvbuf,recvcount,MPI_DOUBLE_PRECISION
     &               ,RANK+block-diag,RANK+block
     &        ,MPI_COMM_WORLD,STATUS,IERROR)


! Place the newly recieved data into a temporary storage array
        DO I=0,NX_S
          DO K=0,NZ_S
            DO J=0,NY_S
              IF (K+(NZ_S+1)*BLOCK.le.NZ) THEN
! Make sure that we don't write beyond the end of the array
                B(I,K+(NZ_S+1)*BLOCK,J)=RECVBUF(I,K,J)
              END IF
            END DO
          END DO
        END DO

        if ((rank.eq.1)) then 
          do i=0,NX_S
          do k=0,NZ_S
            write(*,*) 'block,i,k,B: ',block,i,k,B(i,k,1)
          end do
          end do
        end if

        ELSE
! Here block=diag, we just need to copy locally from A->B

          DO I=0,NX_S
            DO K=0,NZ_S
              DO J=0,NY_S
                sendbuf(I,K,J)=A_X(I+(NX_S+1)*BLOCK,K,J)
              END DO
            END DO
          END DO
          if (rank.eq.1) then
            write(*,*) 'i,k,j,b(0,0): ',i,k,j,b(0,0,1)
          end if
          DO I=0,NX_S
            DO K=0,NZ_S
              DO J=0,NY_S
              IF (K+(NZ_S+1)*BLOCK.le.NZ) THEN
                  B(I,K+(NZ_S+1)*BLOCK,J)=sendbuf(I,K,J)
              END IF
              END DO
            END DO
          END DO
        if ((rank.eq.1)) then 
          do i=0,NX_S
          do k=0,NZ_S
            write(*,*) 'block,i,k,B: ',block,i,k,B(i,k,1)
          end do
          end do
        end if
        END IF

      END DO

! Finally, copy the temporary storage array to the final array
      DO I=0,NX_S
        DO K=0,NZ
          DO J=0,NY_S
            A_Z(I,K,J)=B(I,K,J)
          END DO
        END DO
      END DO

      RETURN
      END       

      SUBROUTINE TRANSPOSE_MPI_X_TO_Y(A_X,A_Y,NX,NY,NZ,NX_S,NY_S,NZ_S)
      INCLUDE "/usr/local/mpich2/include/mpif.h"
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables

      integer NX,NY,NZ,NX_S,NY_S,NZ_S

      complex*16 A_X(0:NX,0:NZ_S,0:NY_S)
      complex*16 A_Y(0:NX_S,0:NZ_S,0:NY)
      complex*16 B(0:NX_S,0:NZ_S,0:NY)

      complex*16 SENDBUF(0:NX_S,0:NZ_S,0:NY_S)
      complex*16 RECVBUF(0:NX_S,0:NZ_S,0:NY_S)

      integer DIAG
      integer block
      integer sendcount,recvcount

      NP_S=NPROCS**0.5

      DIAG=MOD(RANK,NP_S)

      sendcount=(NX_S+1)*(NZ_S+1)*(NY_S+1)
      recvcount=sendcount

      DO block=0,NP_S-1

! Get a block from the source array to be sent to another process
        DO I=0,NX_S
          DO K=0,NZ_S
            DO J=0,NY_S
              SENDBUF(I,K,J)=A_X(I+(NX_S+1)*BLOCK,K,J)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself
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
      INCLUDE "/usr/local/mpich2/include/mpif.h"
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables

      integer NX,NY,NZ,NX_S,NY_S,NZ_S

      complex*16 A_X(0:NX,0:NZ_S,0:NY_S)
      complex*16 A_Y(0:NX_S,0:NZ_S,0:NY)
      complex*16 B(0:NX,0:NZ_S,0:NY_S)

      complex*16 SENDBUF(0:NX_S,0:NZ_S,0:NY_S)
      complex*16 RECVBUF(0:NX_S,0:NZ_S,0:NY_S)

      integer DIAG
      integer block
      integer sendcount,recvcount

      NP_S=NPROCS**0.5

      DIAG=MOD(RANK,NP_S)

      sendcount=(NX_S+1)*(NZ_S+1)*(NY_S+1)
      recvcount=sendcount

      DO block=0,NP_S-1

! Get a block from the source array to be sent to another process
        DO I=0,NX_S
          DO K=0,NZ_S
            DO J=0,NY_S
              SENDBUF(I,K,J)=A_Y(I,K,J+(NY_S+1)*BLOCK)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself
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
              IF (I+(NX_S+1)*BLOCK.le.NX) THEN
! Make sure that we don't write beyond the end of the array
                B(I+(NX_S+1)*BLOCK,K,J)=RECVBUF(I,K,J)
              END IF
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

      SUBROUTINE TRANSPOSE_MPI_Y_TO_Z(A_Y,A_Z,NX,NY,NZ,NX_S,NY_S,NZ_S)
      INCLUDE "/usr/local/mpich2/include/mpif.h"
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables

      integer NX,NY,NZ,NX_S,NY_S,NZ_S

      real*8 A_Z(0:NX_S,0:NZ,0:NY_S)
      real*8 A_Y(0:NX_S,0:NZ_S,0:NY)
      real*8 B(0:NX_S,0:NZ,0:NY_S)

      real*8 SENDBUF(0:NX_S,0:NZ_S,0:NY_S)
      real*8 RECVBUF(0:NX_S,0:NZ_S,0:NY_S)

      integer DIAG
      integer block
      integer sendcount,recvcount

      NP_S=NPROCS**0.5

      DIAG=MOD(RANK,NP_S)

      sendcount=(NX_S+1)*(NZ_S+1)*(NY_S+1)
      recvcount=sendcount

      DO block=0,NP_S-1

! Get a block from the source array to be sent to another process
        DO I=0,NX_S
          DO K=0,NZ_S
            DO J=0,NY_S
              SENDBUF(I,K,J)=A_Y(I,K,J+(NY_S+1)*BLOCK)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself
          CALL MPI_SENDRECV(
     &        sendbuf,sendcount,MPI_DOUBLE_PRECISION
     &               ,RANK+block-diag,RANK+block
     &        ,recvbuf,recvcount,MPI_DOUBLE_PRECISION
     &               ,RANK+block-diag,RANK+block
     &        ,MPI_COMM_WORLD,STATUS,IERROR)

! Place the newly recieved data into a temporary storage array
        DO I=0,NX_S
          DO K=0,NZ_S
            DO J=0,NY_S
              IF (K+(NZ_S+1)*BLOCK.le.NZ) THEN
! Make sure that we don't write beyond the end of the array
                B(I,K+(NZ_S+1)*BLOCK,J)=RECVBUF(I,K,J)
              END IF
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
              IF (K+(NZ_S+1)*BLOCK.le.NZ) THEN
                B(I,K+(NZ_S+1)*BLOCK,J)=sendbuf(I,K,J)
              END IF
              END DO
            END DO
          END DO
        END IF

      END DO

! Finally, copy the temporary storage array to the final array
      DO I=0,NX_S
        DO K=0,NZ
          DO J=0,NY_S
            A_Z(I,K,J)=B(I,K,J)
          END DO
        END DO
      END DO

      RETURN
      END 

      SUBROUTINE TRANSPOSE_MPI_Z_TO_Y(A_Z,A_Y,NX,NY,NZ,NX_S,NY_S,NZ_S)
      INCLUDE "/usr/local/mpich2/include/mpif.h"
      INCLUDE 'header_mpi'
! This subroutine initializes all mpi variables

      integer NX,NY,NZ,NX_S,NY_S,NZ_S

      real*8 A_Z(0:NX_S,0:NZ,0:NY_S)
      real*8 A_Y(0:NX_S,0:NZ_S,0:NY)
      real*8 B(0:NX_S,0:NZ_S,0:NY)

      real*8 SENDBUF(0:NX_S,0:NZ_S,0:NY_S)
      real*8 RECVBUF(0:NX_S,0:NZ_S,0:NY_S)

      integer DIAG
      integer block
      integer sendcount,recvcount

      NP_S=NPROCS**0.5

      DIAG=MOD(RANK,NP_S)

      sendcount=(NX_S+1)*(NZ_S+1)*(NY_S+1)
      recvcount=sendcount

      DO block=0,NP_S-1

! Get a block from the source array to be sent to another process
        DO I=0,NX_S
          DO K=0,NZ_S
            DO J=0,NY_S
              SENDBUF(I,K,J)=A_Z(I,K+(NZ_S+1)*BLOCK,J)
            END DO
          END DO
        END DO

        IF (block.ne.DIAG) then
! Communicate only if we aren't talking to ourself
          CALL MPI_SENDRECV(
     &        sendbuf,sendcount,MPI_DOUBLE_PRECISION
     &               ,RANK+block-diag,RANK+block
     &        ,recvbuf,recvcount,MPI_DOUBLE_PRECISION
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
                sendbuf(I,K,J)=A_Z(I,K+(NZ_S+1)*BLOCK,J)
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


 
