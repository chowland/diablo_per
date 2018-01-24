C******************************************************************************|
C fft.f, the FFT package for diablo.                               VERSION 0.9
C
C This file isolates all calls to the FFTW package (available at: www.fftw.org)
C These wrapper routines were written by T. Bewley (spring 2001).
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C The arrangement of the significant real numbers in the arrays (denoted by +)
C in physical space, in Fourier space, and in Fourier space after packing are
C shown below for the 2D (X-Z) plane.  The third direction (Y) is handled in
C an identical matter as the Z direction shown here.
C
C      oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
C      oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
C NZ-1 ++++++++++++++++oo     -1  ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo     -2  ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo     -3  ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo         ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo    -NKZ ++++++++++++oooooo         oooooooooooooooooo
C      ++++++++++++++++oo         oooooooooooooooooo     -1  ++++++++++++oooooo
C      ++++++++++++++++oo         oooooooooooooooooo     -2  ++++++++++++oooooo
C      ++++++++++++++++oo         oooooooooooooooooo     -3  ++++++++++++oooooo
C      ++++++++++++++++oo         oooooooooooooooooo         ++++++++++++oooooo
C      ++++++++++++++++oo         oooooooooooooooooo    -NKZ ++++++++++++oooooo
C      ++++++++++++++++oo     NKZ ++++++++++++oooooo     NKZ ++++++++++++oooooo
C      ++++++++++++++++oo         ++++++++++++oooooo         ++++++++++++oooooo
C   3  ++++++++++++++++oo      3  ++++++++++++oooooo      3  ++++++++++++oooooo
C   2  ++++++++++++++++oo      2  ++++++++++++oooooo      2  ++++++++++++oooooo
C   1  ++++++++++++++++oo      1  ++++++++++++oooooo      1  ++++++++++++oooooo
C   0  ++++++++++++++++oo      0  +o++++++++++oooooo      0  +o++++++++++oooooo
C      ^^^^           ^           ^ ^ ^     ^                ^ ^ ^     ^
C      0123           NX-1        0 1 2     NKX              0 1 2     NKX
C
C       PHYSICAL SPACE              FOURIER SPACE         FOURIER SPACE (PACKED)
C
C After the Real->Fourier transform, the significant coefficients are put next
C to each other in the array, so a loop such as
C
C        DO K=0,TNKZ           [where TNKZ = 2*NKZ = 2*(NZ/3) ]
C          DO I=0,NKX          [where  NKX = NX/3             ]
C            CP(I,K,J)= ...
C          END DO
C        END DO
C
C includes all the Fourier coefficients of interest.  The subsequent loops in
C Fourier space just work on these coefficients in the matrix.
C  
C Before a Fourier->Real transform, the significant coefficients are unpacked
C and the higher wavenumbers are SET TO ZERO before the inverse transform.
C This has the effect of doing the required dealiasing.
C
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

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


      PI = 4. * ATAN(1.0)
      CI = CMPLX(0.0,1.0)
      EPS= 0.000000001

      IF (NUM_PER_DIR .GT. 0) THEN
        CALL RFFTWND_F77_CREATE_PLAN(FFTW_X_TO_F_PLAN, 1, NX,
     *        FFTW_FORWARD,  FFTW_MEASURE )
        CALL RFFTWND_F77_CREATE_PLAN(FFTW_X_TO_P_PLAN, 1, NX,
     *        FFTW_BACKWARD, FFTW_MEASURE )
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
        CALL RFFTWND_F77_CREATE_PLAN(FFTW_Z_TO_F_PLAN, 1, NZ,
     *       FFTW_FORWARD,  FFTW_MEASURE )
        CALL RFFTWND_F77_CREATE_PLAN(FFTW_Z_TO_P_PLAN, 1, NZ,
     *       FFTW_BACKWARD, FFTW_MEASURE )
        NKZ=NZ/3
        TNKZ=2*NKZ
        RNZ=1.0*NZ
        DO K=0,NKZ
          KZ(K)=-K*(2.*PI)/LZ
        END DO
        DO K=1,NKZ
          KZ(TNKZ+1-K)=K*(2.*PI)/LZ
        END DO
        DO K=0,TNKZ
          KZ2(K)=KZ(K)*KZ(K)
          CIKZ(K)=CI*KZ(K)
        END DO
      END IF

      IF (NUM_PER_DIR .GT. 2) THEN
        CALL RFFTWND_F77_CREATE_PLAN(FFTW_Y_TO_F_PLAN, 1, NY,
     *       FFTW_FORWARD,  FFTW_MEASURE )
        CALL RFFTWND_F77_CREATE_PLAN(FFTW_Y_TO_P_PLAN, 1, NY,
     *       FFTW_BACKWARD, FFTW_MEASURE )
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

      IF (USE_MPI) THEN

        NKX_S=NINT((NKX+1)/NP_S+0.5)-1
        TNKZ_S=NINT((TNKZ+1)/NP_S+0.5)-1
        TNKY_S=NINT((TNKY+1)/NP_S+0.5)-1

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
          IF ((K+INT(RANK/NP_S)*(TNKZ_S+1)).le.TNKZ) THEN
            CIKZ_S(K)=CIKZ(K+INT(RANK/NP_S)*(TNKZ_S+1))
            KZ_S(K)=KZ(K+INT(RANK/NP_S)*(TNKZ_S+1))
            KZ2_S(K)=KZ2(K+INT(RANK/NP_S)*(TNKZ_S+1))
          END IF
        END DO
        DO J=0,TNKY_S
          IF ((J+MOD(RANK,NP_S)*(TNKY_S+1)).le.TNKY) THEN
            CIKY_S(J)=CIKY(J+MOD(RANK,NP_S)*(TNKY_S+1))
            KY_S(J)=KY(J+MOD(RANK,NP_S)*(TNKY_S+1))
            KY2_S(J)=KY2(J+MOD(RANK,NP_S)*(TNKY_S+1))
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

C Intermediate arrays local X, and Z, These are all equivalenced
      COMPLEX*16 CU_X(0:NX/2,0:NZ_S,0:NY_S) 
      REAL*8 U_X(0:NX+1,0:NZ_S,0:NY_S)
      REAL*8 U_Y(0:NX_S,0:NZ_S,0:NY+1)
      REAL*8 U_Z(0:NX_S,0:NZ+1,0:NY_S)
 
C Equivalence the intermediate arrrays to avoid wasting memory
C The FFTs are done in-place, so this is safe 
      EQUIVALENCE(U_X,U_Y,U_Z)

      DO I=0,NX/2
        DO K=0,NZ_S
          DO j=0,NY_S
            CU_X(I,K,J)=CU(I,K,J)
          END DO
        END DO
      END DO

C Inverse transform in the x-direction:
      CALL FFT_X_TO_PHYSICAL(CU_X,U_X)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C In physical space, we are now indexed from 0..NX-1
C Set the rows above this to zero
      do j=0,NY_S
      do k=0,NZ_S
      do i=NX,NX+1 
        U_X(i,k,j)=0.d0
      end do
      end do
      end do

C Perform a parallel transpose to get data stored locally in the z-direction
      CALL TRANSPOSE_MPI_X_TO_Z(U_X,U_Z,NX+1,NY+1,NZ+1,NX_S,NY_S,NZ_S)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C Inverse transform in the z-direction:
      CALL FFT_Z_TO_PHYSICAL(U_Z,U_Z)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C In Physical space, we are now indexed from 0..NZ-1
C Set the rows above this to zero
      do j=0,NY_S
      do k=NZ,NZ+1
      do i=0,NX_S
        U_Z(i,k,j)=0.d0
      end do
      end do
      end do

C Perform a parallel transpose to get data stored locally in the y-direction
      CALL TRANSPOSE_MPI_Z_TO_Y(U_Z,U_Y,NX+1,NY+1,NZ+1,NX_S,NY_S,NZ_S)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C Inverse transform in the y-direction:
      CALL FFT_Y_TO_PHYSICAL(U_Y,U_Y)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C In Physical space, we are now indexed from 0..NY-1
C Set the rows above this to zero
      do j=NY,NY+1
      do k=0,NZ_S
      do i=0,NX_S
        U_Y(i,k,j)=0.d0
      end do
      end do
      end do

C Now, transfer from U_Y to the output array
      DO I=0,NX_S
        DO K=0,NZ_S
          DO J=0,NY+1
            U(I,K,J)=U_Y(I,K,J)
          END DO
        END DO
      END DO

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
      COMPLEX*16 CU(0:NX/2,0:NZ_S,0:NY_S)

C Intermediate arrays local X, and Z
      REAL*8 U_X(0:NX+1,0:NZ_S,0:NY_S)
      COMPLEX*16 CU_X(0:NX/2,0:NZ_S,0:NY_S)
      REAL*8 U_Y(0:NX_S,0:NZ_S,0:NY+1)
      REAL*8 U_Z(0:NX_S,0:NZ+1,0:NY_S)

      REAL*8 RNUM

C Equivalence the intermediate arrrays to avoid wasting memory
C The FFTs are done in-place, so this is safe
      EQUIVALENCE(U_X,U_Y,U_Z)

      CALL RANDOM_SEED


C Perform an FFT in the y-direction:
      CALL FFT_Y_TO_FOURIER(U,U)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)


! Copy U into the dummy array U_Y which is equivalenced to U_Z
       do i=0,NX_S
         do k=0,NZ_S
           do j=0,NY+1
             U_Y(i,k,j)=U(i,k,j)
           end do
        end do
       end do

C Parallel transpose to get data stored locally in the z-direction
      CALL TRANSPOSE_MPI_Y_TO_Z(U_Y,U_Z,NX+1,NY+1,NZ+1,NX_S,NY_S,NZ_S)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      do k=0,NZ+1
        write(710+RANK,*) k,U_Z(1,k,0)
      end do
C FFT in the z-direction:
      CALL FFT_Z_TO_FOURIER(U_Z,U_Z)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      do k=0,NZ+1
        write(720+RANK,*) k,U_Z(1,k,0)
      end do
C Parallel transpose to get data stored locally in the x-direction
      CALL TRANSPOSE_MPI_Z_TO_X(U_Z,U_X,NX+1,NY+1,NZ+1,NX_S,NY_S,NZ_S)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C FFT in the x-direction:

      do i=0,NX+1
        write(725+RANK,*) I,U_X(I,1,0)
      end do

      CALL FFT_X_TO_FOURIER(U_X,CU_X)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

      do i=0,NX/2
        write(775+RANK,*) I,CU_X(I,1,0)
      end do

C Now move the data to the output array
      

      DO I=0,NX/2
        DO K=0,NZ_S
          DO J=0,NY_S
            CU(I,K,J)=CU_X(I,K,J)
          END DO
        END DO
      END DO 

C Multiply the imaginary parts of the y and z arrays by sqrt(-1)
!      DO I=0,NX/2
!        DO K=0,NZ_S
!          DO J=0,NY_S
!            if (j+mod(RANK,NP_S)*(TNKY_S+1).ge.NKY+1) then
!              CU(I,K,J)=CI*CU(I,K,J)
!            end if
!          END DO
!        END DO
!      END DO
!      DO I=0,NX/2
!        DO K=0,NZ_S
!          DO J=0,NY_S
!            if (k+int(RANK/NP_S)*(TNKZ_S+1).ge.NKZ+1) then 
!              CU(I,K,J)=CI*CU(I,K,J)
!            end if
!          END DO
!        END DO
!      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Z_TO_PHYSICAL(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms along the z direction to Fourier space
C The input and output are real arrays with the input array packed
C to contain both the real and imagingary parts of the transformed array
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER  I, J, K
      REAL*8   CU(0:NX_S,0:NZ+1,0:NY_S)
      REAL*8   U (0:NX_S,0:NZ+1,0:NY_S)

      REAL*8 RZX_PLANE(0:NZ+1,0:NX_S)
      COMPLEX*16 CZX_PLANE(0:NZ/2,0:NX_S)

C First, unpack data into the CZX_PLANE complex temporary
C storage variable (dealiasing in the process)
C Then, perform a complex -> real transform in the z-direction

      DO J=0,NY_S
        DO I=0,NX_S
! Zero CZX_PLANE
          DO K=1,NZ+1
            CZX_PLANE(K,I)=DCMPLX(0.0,0.0)
          END DO
! Mean mode contains only a real part
          K=0
          CZX_PLANE(K,I)=DCMPLX(CU(I,K,J),0.d0)
! All nonzero modes:
          DO K=1,NKZ
            CZX_PLANE(K,I)=DCMPLX(CU(I,K,J),CU(I,TNKZ+1-K,J))
          END DO
        END DO

        CALL RFFTWND_F77_COMPLEX_TO_REAL(FFTW_Z_TO_P_PLAN, NX_S+1,
     *    CZX_PLANE(0,0), 1, NZ/2+1, RZX_PLANE(0,0), 1, NZ+2)

        DO I=0,NX_S
          DO K=0,NZM
            U(I,K,J)=RZX_PLANE(K,I)
          END DO
        END DO
       END DO

       RETURN
       END 

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Z_TO_FOURIER(U,CU)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms in the z-direction
C The input and output should be in physical space with the output
C packed to hold the real and imaginary parts
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE '/usr/local/mpich2/include/mpif.h'
      INCLUDE 'header'
      INCLUDE 'header_mpi'
      INTEGER I, J, K
      REAL*8     U (0:NX_S,0:NZ+1,0:NY_S)
      REAL*8     CU(0:NX_S,0:NZ+1,0:NY_S)

      REAL*8 RZX_PLANE(0:NZ+1,0:NX_S)
      COMPLEX*16 CZX_PLANE(0:NZ/2,0:NX_S)

C First, put the data into the CZX_PLANE temporary storage variable,
C perform a real -> complex transform in the z direction.

      DO J=0,NY_S
        DO I=0,NX_S
          DO K=0,NZ+1
            RZX_PLANE(K,I)=U(I,K,J)
          END DO
        END DO

        CALL RFFTWND_F77_REAL_TO_COMPLEX(FFTW_Z_TO_F_PLAN, NX_S+1,
     *    RZX_PLANE(0,0), 1, NZ+2, CZX_PLANE(0,0), 1, NZ/2+1)

C Scale by NZ (necessary by FFTW convention)

        if (J.eq.0) then
          do k=0,NZ+1
            write(825+RANK,*) k,CZX_PLANE(k,1)/RNZ
          end do
        end if

        DO I=0,NX_S
          DO K=0,NZ/2
            CZX_PLANE(K,I)=CZX_PLANE(K,I)/RNZ
          END DO
        END DO

        DO I=0,NX_S
          DO K=0,NKZ
            CU(I,K,J)=DBLE(CZX_PLANE(K,I))
          END DO
          DO K=1,NKZ
            CU(I,TNKZ+1-K,J)=DIMAG(CZX_PLANE(K,I))
          END DO
        END DO
      END DO


      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Y_TO_FOURIER(U,CU)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms in the y-direction
C The input and output should be in physical space with the output
C packed to hold the real and imaginary parts
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE '/usr/local/mpich2/include/mpif.h'
      INCLUDE 'header'
      INCLUDE 'header_mpi'
      INTEGER I, J, K
      REAL*8     U (0:NX_S,0:NZ_S,0:NY+1)
      REAL*8     CU(0:NX_S,0:NZ_S,0:NY+1)

      REAL*8 RYX_PLANE(0:NY+1,0:NX_S)
      COMPLEX*16 CYX_PLANE(0:NY/2,0:NX_S)

C First, put the data into the CYX_PLANE temporary storage variable,
C perform a real -> complex transform in the y direction.

      DO K=0,NZ_S
        DO I=0,NX_S
          DO J=0,NY+1
            RYX_PLANE(J,I)=U(I,K,J)
          END DO
        END DO

        CALL RFFTWND_F77_REAL_TO_COMPLEX(FFTW_Y_TO_F_PLAN, NX_S+1,
     *    RYX_PLANE(0,0), 1, NY+2, CYX_PLANE(0,0), 1, NY/2+1)

C Scale by NY (necessary by FFTW convention)
        DO I=0,NX_S
          DO J=0,NY/2
            CYX_PLANE(J,I)=CYX_PLANE(J,I)/RNY
          END DO
        END DO            

        DO I=0,NX_S
          DO J=0,NKY
            CU(I,K,J)=DBLE(CYX_PLANE(J,I))
          END DO
          DO J=1,NKY
            CU(I,K,TNKY+1-J)=DIMAG(CYX_PLANE(J,I))
          END DO
        END DO
      END DO

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Y_TO_PHYSICAL(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms along the y direction to Fourier space
C The input and output are real arrays with the input array packed
C to contain both the real and imagingary parts of the transformed array
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER  I, J, K
      REAL*8   CU(0:NX_S,0:NZ_S,0:NY+1)
      REAL*8   U (0:NX_S,0:NZ_S,0:NY+1)

      REAL*8 RYX_PLANE(0:NY+1,0:NX_S)
      COMPLEX*16 CYX_PLANE(0:NY/2,0:NX_S)

C First, unpack data into the CYX_PLANE complex temporary
C storage variable (dealiasing in the process)
C Then, perform a complex -> real transform in the y-direction

      DO K=0,NZ_S
        DO I=0,NX_S
! Zero CYX_PLANE
          DO J=1,NY+1
            CYX_PLANE(J,I)=DCMPLX(0.0,0.0)
          END DO
! Mean mode contains only a real part
          J=0
          CYX_PLANE(J,I)=DCMPLX(CU(I,K,J),0.d0)
! All nonzero modes:
          DO J=1,NKY
            CYX_PLANE(J,I)=DCMPLX(CU(I,K,J),CU(I,K,TNKY+1-J))
          END DO
        END DO

        CALL RFFTWND_F77_COMPLEX_TO_REAL(FFTW_Y_TO_P_PLAN, NX_S+1,
     *    CYX_PLANE(0,0), 1, NY/2+1, RYX_PLANE(0,0), 1, NY+2)

        DO I=0,NX_S
          DO J=0,NY+1
            U(I,K,J)=RYX_PLANE(J,I)
          END DO
        END DO
       END DO

       RETURN
       END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_X_TO_FOURIER(U,CU)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I, J, K
      REAL*8     U (0:NX+1,0:NZ_S,0:NY_S)
      COMPLEX*16 CU(0:NX/2,0:NZ_S,0:NY_S)

C Looping over the planes of interest, simply perform a real -> complex
C transform in place in the big storage array, scaling appropriately.

      DO J=0,NY_S
       CALL RFFTWND_F77_REAL_TO_COMPLEX(FFTW_X_TO_F_PLAN,(NZ_S+1),
     *    U(0,0,J), 1, NX+2, CU(0,0,J), 1, NX/2+1)
        DO K=0,NZ_S
          DO I=0,NKX
            CU(I,K,J)=CU(I,K,J)/RNX
          END DO
        END DO
      END DO

      do i=0,NX/2
        write(750+RANK,*) I,CU(I,1,0)
      end do

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_X_TO_PHYSICAL(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) to physical space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE '/usr/local/mpich2/include/mpif.h'
      INCLUDE 'header'
      INCLUDE 'header_mpi'

      INTEGER I, J, K
      REAL*8     U (0:NX+1,0:NZ_S,0:NY_S)
      COMPLEX*16 CU(0:NX/2,0:NZ_S,0:NY_S)

C Looping over the planes of interest, simply set the higher wavenumbers to
C zero and then perform a complex -> real transform in place in the big
C storage array.

      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=NKX+1,NX/2
C Perform De-aliasing on the largest third of the wavenumbers
            CU(I,K,J)=0.
          END DO
        END DO

       CALL RFFTWND_F77_COMPLEX_TO_REAL(FFTW_X_TO_P_PLAN,(NZ_S+1),
     *    CU(0,0,J), 1, NX/2+1, U(0,0,J), 1, NX+2)

      END DO


      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZ_TO_FOURIER(U,CU,JMIN,JMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 2 directions) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZ_TO_PHYSICAL(CU,U,JMIN,JMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 2 directions) planes JMIN-JMAX to physical space.
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


      RETURN
      END





