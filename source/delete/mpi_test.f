

      program mpi_test

      INCLUDE "/usr/local/mpich2/include/mpif.h"
      INCLUDE 'header_mpi'
      INCLUDE 'header'
! This subroutine initializes all mpi variables

!      integer NX,NY,NZ
!      parameter (NX=64,NY=64,NZ=64)
!      parameter (NX_S=16,NY_S=16,NZ_S=16)

      real*8 A(0:NX_S,0:NZ_S,0:NY+1)
      complex*16 B(0:NX/2,0:NZ_S,0:NY_S)

      real*8 A_X(0:NX+1,0:NZ_S,0:NY_S)
      real*8 A_Y(0:NX_S,0:NZ_S,0:NY+1)
      real*8 A_Z(0:NX_s,0:NZ+1,0:NY_S)

      EQUIVALENCE(A,B,A_X,A_Y,A_Z)

      CALL MPI_INIT(IERROR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERROR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,RANK,IERROR)

      CALL INIT_FFT

      do i=0,NX_S
        do k=0,NZ_S
          do j=0,NY
            A(i,k,j)=i+j+k
          end do
        end do
      end do

        k=1
        do i=0,NX_S
        do j=0,NY
          write(100+rank,*) 'i,j,A(i,k): ',i,j,A(i,k,j)
        end do
        end do

      CALL FFT_XZY_MPI_TO_PHYSICAL(A,B)
      CALL FFT_XZY_MPI_TO_FOURIER(B,A)


      k=1
        do i=0,NX
          do j=0,NY_S
            write(200+rank,*) 'i,j,A(i,k): ',i,j,A(i,k,j)
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


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_FFT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE "/usr/local/mpich2/include/mpif.h"

      INCLUDE 'header'
      INCLUDE 'header_mpi'
      INTEGER I,J,K

      INTEGER         FFTW_FORWARD,      FFTW_BACKWARD,
     *                FFTW_ESTIMATE,     FFTW_MEASURE,
     *                FFTW_OUT_OF_PLACE, FFTW_IN_PLACE,
     *                FFTW_USE_WISDOM,   FFTW_THREADSAFE
      PARAMETER(      FFTW_FORWARD=-1,      FFTW_BACKWARD=1,
     *                FFTW_ESTIMATE=0,      FFTW_MEASURE=1,
     *                FFTW_OUT_OF_PLACE=0,  FFTW_IN_PLACE=8,
     *                FFTW_USE_WISDOM=16,   FFTW_THREADSAFE=128 )

      WRITE(6,*) 'Initializing FFTW package.'
!      write(*,*) 'In init FFT, RANK,NPROCS: ',RANK,NPROCS

      write(*,*) 'done writing test'

      PI = 4. * ATAN(1.0)
      CI = CMPLX(0.0,1.0)
      EPS= 0.000000001

      IF (NUM_PER_DIR .GT. 0) THEN
        CALL RFFTWND_F77_CREATE_PLAN(FFTW_X_TO_F_PLAN, 1, NX,
     *        FFTW_FORWARD,  FFTW_MEASURE + FFTW_IN_PLACE)
        CALL RFFTWND_F77_CREATE_PLAN(FFTW_X_TO_P_PLAN, 1, NX,
     *        FFTW_BACKWARD, FFTW_MEASURE + FFTW_IN_PLACE)
        CALL RFFTW_F77_CREATE_PLAN(FFTW_TEST_PLAN,NX,FFTW_BACKWARD
     &           ,FFTW_MEASURE )
        NKX=NX/3
        RNX=1.0*NX
        DO I=0,NKX
          KX(I)=I*(2.*PI)/LX
          KX2(I)=KX(I)*KX(I)
          CIKX(I)=CI*KX(I)
        END DO
 
      END IF

      IF (NUM_PER_DIR .GT. 1) THEN
        CALL FFTWND_F77_CREATE_PLAN(FFTW_Z_TO_F_PLAN, 1, NZ,
     *       FFTW_FORWARD,  FFTW_MEASURE  + FFTW_IN_PLACE)
        CALL FFTWND_F77_CREATE_PLAN(FFTW_Z_TO_P_PLAN, 1, NZ,
     *       FFTW_BACKWARD, FFTW_MEASURE  + FFTW_IN_PLACE)
        NKZ=NZ/3
        TNKZ=2*NKZ
        RNZ=1.0*NZ
        DO K=0,NKZ
          KZ(K)=K*(2.*PI)/LZ
        END DO
        DO K=1,NKZ
          KZ(TNKZ+1-K)=-K*(2.*PI)/LZ
        END DO
        DO K=0,TNKZ
          KZ2(K)=KZ(K)*KZ(K)
          CIKZ(K)=CI*KZ(K)
        END DO
      END IF

      IF (NUM_PER_DIR .GT. 2) THEN
        CALL FFTWND_F77_CREATE_PLAN(FFTW_Y_TO_F_PLAN, 1, NY,
     *       FFTW_FORWARD,  FFTW_MEASURE + FFTW_IN_PLACE)
        CALL FFTWND_F77_CREATE_PLAN(FFTW_Y_TO_P_PLAN, 1, NY,
     *       FFTW_BACKWARD, FFTW_MEASURE + FFTW_IN_PLACE)
        NKY=NY/3
        TNKY=2*NKY
        RNY=1.0*NY
	  KY(0) = 0.
	  KY2(0) = 0.
	  CIKY(0) = (0.0,0.0)
        DO J=1,NKY
          KY(J)=J*(2.*PI)/LY
        END DO
        DO J=1,NKY
          KY(TNKY+1-J)=-J*(2.*PI)/LY
        END DO
        DO J=1,TNKY
          KY2(J)=KY(J)*KY(J)
          CIKY(J)=CI*KY(J)
        END DO
      END IF

      DO I=0,NKX
        DO K=0,NZ
          CZX_PLANE(K,I)=CMPLX(0.0,0.0)
        END DO
      END DO
      DO K=0,TNKZ
        DO J=0,NY
          CYZ_PLANE(J,K)=CMPLX(0.0,0.0)
        END DO
      END DO


      IF (USE_MPI) THEN

        WRITE(*,*) 'Checkpoint',RANK,NPROCS

        NKX_S=NINT((NKX+1)/NPROCS+0.5)-1
        TNKZ_S=NINT((TNKZ+1)/NPROCS+0.5)-1
        TNKY_S=NINT((TNKY+1)/NPROCS+0.5)-1

        WRITE(*,*) 'NKX,TNKZ,TNKY',NKX,TNKZ,TNKY
        write(*,*) 'NPROCS: ',NPROCS
        WRITE(*,*) 'NKX_S,TNKZ_S,TNKY_S',NKX_S,TNKZ_S,TNKY_S
        write(*,*) 'NX_S,NY_S,NZ_S: ',NX_S,NY_S,NZ_S

        DO I=0,NKX_S
          IF ((I+RANK*(NKX_S+1)).le.NKX) THEN
            CIKX_S(I)=CIKX(I+RANK*(NKX_S+1))
            KX_S(I)=KX(I+RANK*(NKX_S+1))
            KX2_S(I)=KX2(I+RANK*(NKX_S+1))
          END IF
        END DO
        DO K=0,TNKZ_S
          IF ((K+RANK*(TNKZ_S+1)).le.TNKZ) THEN
            CIKZ_S(K)=CIKZ(K+RANK*(TNKZ_S+1))
            KZ_S(K)=KZ(K+RANK*(TNKZ_S+1))
            KZ2_S(K)=KZ2(K+RANK*(TNKZ_S+1))
          END IF
        END DO
        DO J=0,TNKY_S
          IF ((J+RANK*(TNKY_S+1)).le.TNKY) THEN
            CIKY_S(J)=CIKY(J+RANK*(TNKY_S+1))
            KY_S(J)=KY(J+RANK*(TNKY_S+1))
            KY2_S(J)=KY2(J+RANK*(TNKY_S+1))
          END IF
        END DO
      END IF


      WRITE(6,*) 'FFTW package initialized.'

      RETURN
      END

C******************************************************************************|
C-------------> The transform routines for the duct flow follow. <-------------|
C******************************************************************************|


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZY_MPI_TO_PHYSICAL(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      include "/usr/local/mpich2/include/mpif.h"
      include 'header_mpi'

      INTEGER JMIN, JMAX, KMIN, KMAX, I, J, K
C Input array is in Fourier space and local in X
      COMPLEX*16 CU(0:NX/2,0:NZ_S,0:NY_S)
C Output array is in physical space and local in Y
      REAL*8     U (0:NX_S,0:NZ_S,0:NY+1)

C Intermediate arrays local X, and Z
      REAL*8 U_X(0:NX+1,0:NZ_S,0:NY_S)
      REAL*8 U_Y(0:NX_S,0:NZ_S,0:NY+1)
      REAL*8 U_Z(0:NX_S,0:NZ+1,0:NY_S)
 
C Equivalence the intermediate arrrays to avoid wasting memory
C The FFTs are done in-place, so this is safe 
      EQUIVALENCE(U_X,U_Y,U_Z)


C Inverse transform in the x-direction:
      CALL FFT_X_TO_PHYSICAL(CU,U_X,0,NY_S,0,NZ_S)
      write(*,*) 'IN FFT: RANK,U_X: ',RANK,U_X(5,5,5)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C Perform a parallel transpose to get data stored locally in the z-direction
      CALL TRANSPOSE_MPI_X_TO_Z(U_X,U_Z,NX,NY,NZ,NX_S,NY_S,NZ_S)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C Inverse transform in the z-direction:
      CALL FFT_Z_TO_PHYSICAL(U_Z,U_Z,0,NY_S)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C Perform a parallel transpose to get data stored locally in the y-direction
      CALL TRANSPOSE_MPI_Z_TO_Y(U_Z,U_Y,NX,NY,NZ,NX_S,NY_S,NZ_S)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C Inverse transform in the y-direction:
      CALL FFT_Y_TO_PHYSICAL(U_Y,U_Y)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)


C Now, transfer from U_Y to the output array
      DO I=0,NX_S
        DO K=0,NZ_S
          DO J=0,NY+1
            U(I,K,J)=U_Y(I,K,J)
          END DO
        END DO
      END DO

      write(*,*) 'DONE FFT'

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZY_MPI_TO_FOURIER(U,CU)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      include "/usr/local/mpich2/include/mpif.h"
      include 'header_mpi'
      
      INTEGER JMIN, JMAX, KMIN, KMAX, I, J, K
C Input array is in Physical space and local in Y
      REAL*8     U (0:NX_S,0:NZ_S,0:NY+1)
C Output array is in Fourier space and local in X
      COMPLEX*16 CU(0:NX+1,0:NZ_S,0:NY_S)

C Intermediate arrays local X, and Z
      REAL*8 U_X(0:NX+1,0:NZ_S,0:NY_S)
      REAL*8 U_Y(0:NX_S,0:NZ_S,0:NY+1)
      REAL*8 U_Z(0:NX_S,0:NZ+1,0:NY_S)

C Equivalence the intermediate arrrays to avoid wasting memory
C The FFTs are done in-place, so this is safe
      EQUIVALENCE(U_X,U_Y,U_Z)

C Perform an FFT in the y-direction:
      CALL FFT_Y_TO_FOURIER(U,U)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C Parallel transpose to get data stored locally in the z-direction
      CALL TRANSPOSE_MPI_Y_TO_Z(U,U_Z,NX,NY,NZ,NX_S,NY_S,NZ_S)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C FFT in the z-direction:
      CALL FFT_Z_TO_FOURIER(U_Z,U_Z,0,NY_S)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C Parallel transpose to get data stored locally in the x-direction
      CALL TRANSPOSE_MPI_Z_TO_X(U_Z,U_X,NX,NY,NZ,NX_S,NY_S,NZ_S)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C FFT in the x-direction:
      CALL FFT_X_TO_FOURIER(U_X,U_X,0,NY_S,0,NZ_S)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C Now, move the data to the output array
      DO I=0,NX+1
        DO K=0,NZ_S
          DO J=0,NY_S
            CU(I,K,J)=U_X(I,K,J)
          END DO
        END DO
      END DO 

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Z_TO_PHYSICAL(CU,U,JMIN,JMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER IMIN, IMAX, JMIN, JMAX, I, J, K
      REAL*8     CU(0:NX_S,0:NZ+1,0:NY_S)
      REAL*8     U (0:NX_S,0:NZ+1,0:NY_S)
    
C Here, start with an array that is physical in x and fourier in z and y
C The data is stored in a real array with the complex fourier amplitudes
C in y and z packed in the way described at the top of fft.f.
C Looping over the planes of interest, unpack data into the CZX_PLANE temporary
C storage variable (setting higher wavenumbers to zero), perform a complex ->
C real transform in the z direction

      DO J=JMIN,JMAX
       IF (NZ.GT.1) THEN
        DO I=0,NX_S
          DO K=0,NKZ
            CZX_PLANE(K,I)=CU(I,K,J)
          END DO
          DO K=NKZ+1,NZM-NKZ
            CZX_PLANE(K,I)=CMPLX(0.0,0.0)
          END DO
          DO K=1,NKZ
            CZX_PLANE(NZM-NKZ+K,I)=CU(I,NKZ+K,J)
          END DO
        END DO
        CALL FFTWND_F77(FFTW_Z_TO_P_PLAN, NX_S,
     *    CZX_PLANE(0,0), 1, NZ+1, CZX_PLANE(0,0), 1, NZ+1)
        DO I=0,NX_S
          DO K=0,NZM
            U(I,K,J)=real(CZX_PLANE(K,I))
          END DO
        END DO
       END IF
       END DO

       RETURN
       END 

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Z_TO_FOURIER(U,CU,JMIN,JMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 2 directions) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER JMIN, JMAX, I, J, K
      REAL*8     U (0:NX_S,0:NZ_S,0:NY+1)
      REAL*8     CU(0:NX_S,0:NZ_S,0:NY+1)

C Looping over the planes of interest, transform in the X direction using
C FFT_X_TO_FOURIER, put the data into the CZX_PLANE temporary storage variable,
C perform a complex -> complex transform in the z direction, then put the data
C back into the big storage array, packing the data towards K=0 and scaling
C appropriately.

      DO J=JMIN,JMAX
        DO I=0,NX_S
          DO K=0,NZM
            CZX_PLANE(J,I)=U(I,K,J)
          END DO
        END DO
        CALL FFTWND_F77(FFTW_Z_TO_F_PLAN, NKX+1,
     *    CZX_PLANE(0,0), 1, NZ+1, CZX_PLANE(0,0), 1, NZ+1)
        DO I=0,NX_S
          DO K=0,NKZ
            CU(I,K,J)=CZX_PLANE(K,I)/RNZ
          END DO
          DO K=1,NKZ
            CU(I,NKZ+K,J)=CZX_PLANE(NZM-NKZ+K,I)/RNZ
          END DO
        END DO
      END DO

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Y_TO_FOURIER(U,CU)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 2 directions) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER JMIN, JMAX, I, J, K
      REAL*8     U (0:NX_S,0:NZ_S,0:NY+1)
      REAL*8     CU(0:NX_S,0:NZ_S,0:NY+1)

C Looping over the planes of interest, transform in the X direction using
C FFT_X_TO_FOURIER, put the data into the CZX_PLANE temporary storage variable,
C perform a complex -> complex transform in the z direction, then put the data
C back into the big storage array, packing the data towards K=0 and scaling
C appropriately.

      DO K=0,NZ_S
        DO I=0,NX_S
          DO J=0,NYM
            CYX_PLANE(J,I)=U(I,K,J)
          END DO
        END DO
        CALL FFTWND_F77(FFTW_Y_TO_F_PLAN, NKX+1,
     *    CYX_PLANE(0,0), 1, NY+1, CYX_PLANE(0,0), 1, NY+1)
        DO I=0,NX_S
          DO J=0,NKY
            CU(I,K,J)=CYX_PLANE(J,I)/RNY
          END DO
          DO J=1,NKY
            CU(I,K,NKY+J)=CYX_PLANE(NYM-NKY+J,I)/RNY
          END DO
        END DO
      END DO

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Y_TO_PHYSICAL(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER IMIN, IMAX, JMIN, JMAX, I, J, K
      REAL*8     CU(0:NX_S,0:NZ_S,0:NY+1)
      REAL*8     U (0:NX_S,0:NZ_S,0:NY+1)

C Here, start with an array that is physical in x and fourier in z and y
C The data is stored in a real array with the complex fourier amplitudes
C in y and z packed in the way described at the top of fft.f.
C Looping over the planes of interest, unpack data into the CZX_PLANE temporary
C storage variable (setting higher wavenumbers to zero), perform a complex ->
C real transform in the z direction

      DO K=0,NZ_S
        DO I=0,NX_S
          DO J=0,NKY
            CYX_PLANE(J,I)=CU(I,K,J)
          END DO
          DO J=NKY+1,NYM-NKY
            CYX_PLANE(J,I)=CMPLX(0.0,0.0)
          END DO
          DO J=1,NKY
            CYX_PLANE(NYM-NKY+J,I)=CU(I,K,NKY+J)
          END DO
        END DO
        CALL FFTWND_F77(FFTW_Y_TO_P_PLAN, NX_S,
     *    CYX_PLANE(0,0), 1, NY+1, CZX_PLANE(0,0), 1, NY+1)
        DO I=0,NX_S
          DO J=0,NYM
            U(I,K,J)=real(CZX_PLANE(J,I))
          END DO
        END DO
       END DO

       RETURN
       END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_X_TO_FOURIER(U,CU,JMIN,JMAX,KMIN,KMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER JMIN, JMAX, KMIN, KMAX, I, J, K
      REAL*8     U (0:NX+1,0:NZ_S,0:NY_S)
      COMPLEX*16 CU(0:NX/2,0:NZ_S,0:NY_S)

C Looping over the planes of interest, simply perform a real -> complex
C transform in place in the big storage array, scaling appropriately.

      DO J=JMIN,JMAX
       CALL RFFTWND_F77_REAL_TO_COMPLEX(FFTW_X_TO_F_PLAN,(KMAX-KMIN+1),
     *    U(0,KMIN,J), 1, NX+2, CU(0,KMIN,J), 1, NX/2+1)
        DO K=KMIN,KMAX
          DO I=0,NKX
            CU(I,K,J)=CU(I,K,J)/RNX
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_X_TO_PHYSICAL(CU,U,JMIN,JMAX,KMIN,KMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes JMIN-JMAX to physical space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE '/usr/local/mpich2/include/mpif.h'
      INCLUDE 'header'
      INCLUDE 'header_mpi'

      INTEGER JMIN, JMAX, KMIN, KMAX, I, J, K
      REAL*8     U (0:NX+1,0:NZ_S,0:NY_S)
      COMPLEX*16 CU(0:NX/2,0:NZ_S,0:NY_S)

      real*8 test(0:NX+1)
      complex*16 ctest(0:NX/2)
!      equivalence (test,ctest)

C Looping over the planes of interest, simply set the higher wavenumbers to
C zero and then perform a complex -> real transform in place in the big
C storage array.

      write(*,*) 'IN FFT_X: NKX: ',NKX
      write(*,*) 'IN FFT_X: CU: ',CU(5,5,5)
      write(*,*) 'IN FFT_X, JMIN,JMAX: ',JMIN,JMAX
      write(*,*) 'IN FFT_X, KMIN,KMAX: ',KMIN,KMAX
      write(*,*) 'IN FFT_X, NY_S,NZ_S: ',NY_S,NZ_S

      write(*,*)  'FFTW_X_TO_P_PLAN: ',FFTW_X_TO_P_PLAN
      write(*,*)  'FFTW_TEST_PLAN: ',FFTW_TEST_PLAN

       do I=0,NX/2
         CTEST(I)=CU(I,5,5)
         write(200+MPI_NUM,*) I,CTEST(I)
         CALL RFFTW_F77_ONE(FFTW_TEST_PLAN,CTEST,TEST)
         write(300+MPI_NUM,*) I,TEST(I),CTEST(I)
       END DO

      DO J=JMIN,JMAX
        DO K=KMIN,KMAX
          DO I=NKX+1,NX/2
C Perform De-aliasing on the largest third of the wavenumbers
            CU(I,K,J)=0.
          END DO
        END DO

       CALL RFFTWND_F77_COMPLEX_TO_REAL(FFTW_X_TO_P_PLAN,(KMAX-KMIN+1),
     *    CU(0,KMIN,J), 1, NX/2+1, U(0,KMIN,J), 1, NX+2)

      END DO

       write(*,*) 'Done FFT_X: U: ',U(5,5,5)

      RETURN
      END

C******************************************************************************|
C-----------> The transform routines for the channel flow follow. <------------|
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZ_TO_FOURIER(U,CU,JMIN,JMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 2 directions) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER JMIN, JMAX, I, J, K
      REAL*8     U (0:NX+1,0:NZ_S,0:NY_S+1)
      COMPLEX*16 CU(0:NX/2,0:NZ_S,0:NY_S+1)

C Looping over the planes of interest, transform in the X direction using
C FFT_X_TO_FOURIER, put the data into the CZX_PLANE temporary storage variable,
C perform a complex -> complex transform in the z direction, then put the data
C back into the big storage array, packing the data towards K=0 and scaling
C appropriately.

      DO J=JMIN,JMAX
        CALL FFT_X_TO_FOURIER(U,CU,J,J,0,NZ-1)
        IF (NZ.GT.1) THEN
        DO I=0,NKX
          DO K=0,NZ-1
            CZX_PLANE(K,I)=CU(I,K,J)
          END DO
        END DO        
        CALL FFTWND_F77(FFTW_Z_TO_F_PLAN, NKX+1,
     *    CZX_PLANE(0,0), 1, NZ+1, CZX_PLANE(0,0), 1, NZ+1)
        DO I=0,NKX
          DO K=0,NKZ
            CU(I,K,J)=CZX_PLANE(K,I)/RNZ
          END DO
          DO K=1,NKZ
            CU(I,NKZ+K,J)=CZX_PLANE(NZ-1-NKZ+K,I)/RNZ
          END DO
        END DO
       END IF
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZ_TO_PHYSICAL(CU,U,JMIN,JMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 2 directions) planes JMIN-JMAX to physical space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER JMIN, JMAX, I, J, K
      REAL*8     U (0:NX+1,0:NZ+1,0:NY+1)
      COMPLEX*16 CU(0:NX/2,0:NZ+1,0:NY+1)

C Looping over the planes of interest, unpack data into the CZX_PLANE temporary
C storage variable (setting higher wavenumbers to zero), perform a complex ->
C complex transform in the z direction, then put the data back into the big
C storage array and transform in the X direction using FFT_X_TO_PHYSICAL.


      DO J=JMIN,JMAX
       IF (NZ.GT.1) THEN
        DO I=0,NKX
          DO K=0,NKZ
            CZX_PLANE(K,I)=CU(I,K,J)
          END DO
          DO K=NKZ+1,NZ-1-NKZ
            CZX_PLANE(K,I)=CMPLX(0.0,0.0)
          END DO
          DO K=1,NKZ
            CZX_PLANE(NZ-1-NKZ+K,I)=CU(I,NKZ+K,J)
          END DO
        END DO
        CALL FFTWND_F77(FFTW_Z_TO_P_PLAN, NKX+1,
     *    CZX_PLANE(0,0), 1, NZ+1, CZX_PLANE(0,0), 1, NZ+1)
        DO I=0,NKX
          DO K=0,NZ-1
            CU(I,K,J)=CZX_PLANE(K,I)
          END DO
        END DO 
       END IF       
       CALL FFT_X_TO_PHYSICAL(CU,U,J,J,0,NZ-1)
      END DO


      RETURN
      END

C******************************************************************************|
C--------> The transform routines for the fully-periodic box follow. <---------|
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZY_TO_FOURIER(U,CU)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 3 directions) the entire box to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I, J, K
      REAL*8     U (0:NX+1,0:NZ+1,0:NY+1)
      COMPLEX*16 CU(0:NX/2,0:NZ+1,0:NY+1)

C First, transform in the X & Z directions using FFT_XZ_TO_FOURIER.  Then,
C looping over the planes of interest, put the data into the CYZ_PLANE
C temporary storage variable, perform a complex -> complex transform in the
C y direction, then put the data back into the big storage array, packing the
C data towards J=0 and scaling appropriately.

      CALL FFT_XZ_TO_FOURIER(U,CU,0,NYM)
      DO I=0,NKX
        DO K=0,TNKZ
          DO J=0,NYM       
            CYZ_PLANE(J,K)=CU(I,K,J)
          END DO
        END DO        
        CALL FFTWND_F77(FFTW_Y_TO_F_PLAN, TNKZ+1,
     *    CYZ_PLANE(0,0), 1, NY+1, CYZ_PLANE(0,0), 1, NY+1)
        DO K=0,TNKZ 
          DO J=0,NKY
            CU(I,K,J)=CYZ_PLANE(J,K)/RNY
          END DO
          DO J=1,NKY
            CU(I,K,NKY+J)=CYZ_PLANE(NYM-NKY+J,K)/RNY
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZY_TO_PHYSICAL(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 3 directions) the entire box to physical space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I, J, K
      REAL*8     U (0:NX+1,0:NZ+1,0:NY+1)
      COMPLEX*16 CU(0:NX/2,0:NZ+1,0:NY+1)

C Looping over the planes of interest, unpack data into the CYZ_PLANE temporary
C storage variable (setting higher wavenumbers to zero), perform a complex ->
C complex transform in the y direction, then put the data back into the big
C storage array.  Finally, transform in the X & Z directions using
C FFT_XZ_TO_PHYSICAL.

      DO I=0,NKX
        DO K=0,TNKZ 
          DO J=0,NKY
            CYZ_PLANE(J,K)=CU(I,K,J)
          END DO
          DO J=NKY+1,NYM-NKY
            CYZ_PLANE(J,K)=CMPLX(0.0,0.0)
          END DO
          DO J=1,NKY
            CYZ_PLANE(NYM-NKY+J,K)=CU(I,K,NKY+J)
          END DO
        END DO
        CALL FFTWND_F77(FFTW_Y_TO_P_PLAN, TNKZ+1,
     *    CYZ_PLANE(0,0), 1, NY+1, CYZ_PLANE(0,0), 1, NY+1)
        DO K=0,TNKZ
          DO J=0,NYM       
            CU(I,K,J)=CYZ_PLANE(J,K)
          END DO
        END DO        
      END DO
      CALL FFT_XZ_TO_PHYSICAL(CU,U,0,NYM)

      RETURN
      END







 
