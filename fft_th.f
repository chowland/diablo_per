C******************************************************************************|
C fft_th.f, the FFT package for diablo.                               VERSION 0.9
C
C This file isolates all calls to the FFTW package (available at: www.fftw.org)
C These wrapper routines were written by T. Bewley (spring 2001).
C This subroutine has been modified to operate on the scalar grid only
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
      SUBROUTINE INIT_FFT_TH
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE "mpif.h"

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

      IF (RANK.eq.0) then
        WRITE(6,*) 'Initializing FFTW package.'
      END IF

      PI = 4. * ATAN(1.0)
      CI = CMPLX(0.0,1.0)
      EPS= 0.000000001

      IF (NUM_PER_DIR .GT. 0) THEN
       CALL RFFTWND_F77_CREATE_PLAN(FFTW_X_TO_F_PLAN_TH, 1, NX_TH,
     *        FFTW_FORWARD,  FFTW_MEASURE )
       CALL RFFTWND_F77_CREATE_PLAN(FFTW_X_TO_P_PLAN_TH, 1, NX_TH,
     *        FFTW_BACKWARD, FFTW_MEASURE )
       CALL RFFTW_F77_CREATE_PLAN(FFTW_TEST_PLAN_TH,NX_TH,FFTW_BACKWARD
     &           ,FFTW_MEASURE )
        NKX_TH=NX_TH/3
        RNX_TH=1.0*NX_TH
        DO I=0,NKX_TH
          KX_TH(I)=I*(2.*PI)/LX
          KX2_TH(I)=KX_TH(I)*KX_TH(I)
          CIKX_TH(I)=CI*KX_TH(I)
        END DO

      END IF

      IF (NUM_PER_DIR .GT. 1) THEN
        CALL FFTWND_F77_CREATE_PLAN(FFTW_Z_TO_F_PLAN_TH, 1, NZ_TH,
     *       FFTW_FORWARD,  FFTW_MEASURE + FFTW_IN_PLACE)
        CALL FFTWND_F77_CREATE_PLAN(FFTW_Z_TO_P_PLAN_TH, 1, NZ_TH,
     *       FFTW_BACKWARD, FFTW_MEASURE + FFTW_IN_PLACE)
        NKZ_TH=NZ_TH/3
        TNKZ_TH=2*NKZ_TH
        RNZ_TH=1.0*NZ_TH
        DO K=0,NKZ_TH
          KZ_TH(K)=K*(2.*PI)/LZ
        END DO
        DO K=1,NKZ_TH
          KZ_TH(TNKZ_TH+1-K)=-K*(2.*PI)/LZ
        END DO
        DO K=0,TNKZ_TH
          KZ2_TH(K)=KZ_TH(K)*KZ_TH(K)
          CIKZ_TH(K)=CI*KZ_TH(K)
        END DO
      END IF

      IF (NUM_PER_DIR .GT. 2) THEN
        CALL FFTWND_F77_CREATE_PLAN(FFTW_Y_TO_F_PLAN_TH, 1, NY_TH,
     *       FFTW_FORWARD,  FFTW_MEASURE + FFTW_IN_PLACE)
        CALL FFTWND_F77_CREATE_PLAN(FFTW_Y_TO_P_PLAN_TH, 1, NY_TH,
     *       FFTW_BACKWARD, FFTW_MEASURE + FFTW_IN_PLACE)
        NKY_TH=NY_TH/3
        TNKY_TH=2*NKY_TH
        RNY_TH=1.0*NY_TH
          KY_TH(0) = 0.
          KY2_TH(0) = 0.
          CIKY_TH(0) = (0.0,0.0)
        DO J=1,NKY_TH
          KY_TH(J)=J*(2.*PI)/LY
        END DO
        DO J=1,NKY_TH
          KY_TH(TNKY_TH+1-J)=-J*(2.*PI)/LY
        END DO
        DO J=1,TNKY_TH
          KY2_TH(J)=KY_TH(J)*KY_TH(J)
          CIKY_TH(J)=CI*KY_TH(J)
        END DO
      END IF

      IF (USE_MPI) THEN

        NKX_S_TH=NINT((NKX_TH+1)/NP_S+0.5)-1
        TNKZ_S_TH=NINT((TNKZ_TH+1)/NP_S+0.5)-1

        DO I=0,NKX_S_TH
          IF ((I+MOD(RANK,NP_S)*(NKX_S_TH+1)).le.NKX_TH) THEN
            CIKX_S_TH(I)=CIKX_TH(I+MOD(RANK,NP_S)*(NKX_S_TH+1))
            KX_S_TH(I)=KX_TH(I+MOD(RANK,NP_S)*(NKX_S_TH+1))
            KX2_S_TH(I)=KX2_TH(I+MOD(RANK,NP_S)*(NKX_S_TH+1))
          END IF
        END DO
        DO K=0,TNKZ_S_TH
          IF ((K+INT(RANK/NP_S)*(TNKZ_S_TH+1)).le.TNKZ_TH) THEN
            CIKZ_S_TH(K)=CIKZ_TH(K+INT(RANK/NP_S)*(TNKZ_S_TH+1))
            KZ_S_TH(K)=KZ_TH(K+INT(RANK/NP_S)*(TNKZ_S_TH+1))
            KZ2_S_TH(K)=KZ2_TH(K+INT(RANK/NP_S)*(TNKZ_S_TH+1))
          END IF
        END DO
      END IF

      IF (RANK.eq.0) THEN
        WRITE(6,*) 'FFTW package initialized.'
      END IF

      RETURN
      END

C******************************************************************************|
C-------------> The transform routines for the duct flow follow. <-------------|
C******************************************************************************|


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZY_MPI_TO_PHYSICAL_TH(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'

      INTEGER JMIN, JMAX, KMIN, KMAX, I, J, K
C Input array is in Fourier space and local in X
      COMPLEX*16 CU(0:NX_S_TH/2,0:NZ_S_TH,0:NY_TH+1)
C Output array is in physical space and local in Y
      REAL*8     U (0:NX_TH+1,0:NZ_S_TH,0:NY_S_TH+1)

C Intermediate arrays local X, and Z, These are all equivalenced
      REAL*8 U_X(0:NX_TH+1,0:NZ_S_TH,0:NY_S_TH)
      COMPLEX*16 CU_X(0:NX_TH/2,0:NZ_S_TH,0:NY_S_TH) 

C CU_Y and CU_Z can't be equivalenced since they are used
C in the transpose routines
      COMPLEX*16 CU_Y(0:NX_S_TH/2,0:NZ_S_TH,0:NY_TH+1)
      COMPLEX*16 CU_Z(0:NX_S_TH/2,0:NZ_TH+1,0:NY_S_TH)
 
C Equivalence the intermediate arrrays to avoid wasting memory
C The FFTs are done in-place, so this is safe 
      EQUIVALENCE(CU_X,U_X)

C Inverse transform in the x-direction:
      CALL FFT_Y_TO_PHYSICAL_TH(CU,CU_Y)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C Perform a parallel transpose to get data stored locally in the z-direction
      CALL TRANSPOSE_MPI_Y_TO_Z_TH(CU_Y,CU_Z)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C Inverse transform in the z-direction:
      CALL FFT_Z_TO_PHYSICAL_TH(CU_Z,CU_Z)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C Perform a parallel transpose to get data stored locally in the y-direction
      CALL TRANSPOSE_MPI_Z_TO_X_TH(CU_Z,CU_X)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C Inverse transform in the y-direction:
      CALL FFT_X_TO_PHYSICAL_TH(CU_X,U_X)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C Now, transfer from U_X to the output array
      DO J=0,NY_S_TH
        DO K=0,NZ_S_TH
          DO I=0,NX_TH+1
            U(I,K,J)=U_X(I,K,J)
          END DO
        END DO
      END DO

C Set data above NXM=0
      DO I=NX_TH,NX_TH+1
        U(I,0:NZ_S,0:NY_S)=0.d0
      END DO
C Set data above NYM=0
      DO J=max(0,NY_TH-INT(RANK/NP_S)*(NY_S_TH+1)),NY_S_TH+1
        U(0:NX_TH+1,0:NZ_S_TH,J)=0.d0
      END DO
C Set data above NZM=0
      DO K=max(0,NZ_TH-MOD(RANK,NP_S)*(NZ_S_TH+1)),NZ_S_TH
        U(0:NX_TH+1,K,0:NY_S_TH)=0.d0
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZY_MPI_TO_FOURIER_TH(U,CU)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'
 
      INTEGER JMIN, JMAX, KMIN, KMAX, I, J, K
C Input array is in Physical space and local in X
      REAL*8     U (0:NX_TH+1,0:NZ_S_TH,0:NY_S_TH+1)
C Output array is in Fourier space and local in Y
      COMPLEX*16 CU(0:NX_S_TH/2,0:NZ_S_TH,0:NY_TH+1)

C Intermediate arrays local X, and Z
      REAL*8 U_X(0:NX_TH+1,0:NZ_S_TH,0:NY_S_TH)
      COMPLEX*16 CU_X(0:NX_TH/2,0:NZ_S_TH,0:NY_S_TH)

      COMPLEX*16 CU_Z(0:NKX_S_TH,0:NZ_TH+1,0:NY_S_TH)
      COMPLEX*16 CU_Y(0:NKX_S_TH,0:TNKZ_S_TH,0:NY_TH+1)

      REAL*8 RNUM

      EQUIVALENCE (U_X,CU_X)

      do j=0,NY_S_TH
      do k=0,NZ_S_TH
      do i=0,NX_TH+1
        U_X(i,k,j)=U(i,k,j)
      end do
      end do
      end do

      CALL FFT_X_TO_FOURIER_TH(U_X,CU_X)

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C Parallel transpose to get data stored locally in the z-direction
      CALL TRANSPOSE_MPI_X_TO_Z_TH(CU_X,CU_Z)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C FFT in the z-direction:
      CALL FFT_Z_TO_FOURIER_TH(CU_Z,CU_Z)

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C Parallel transpose to get data stored locally in the z-direction
      CALL TRANSPOSE_MPI_Z_TO_Y_TH(CU_Z,CU_Y)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C FFT in the y-direction:
      CALL FFT_Y_TO_FOURIER_TH(CU_Y,CU_Y)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

      DO J=0,NY_TH+1
        DO K=0,TNKZ_S_TH
          DO I=0,NKX_S_TH
            CU(I,K,J)=CU_Y(I,K,J)
          END DO
          DO I=NKX_S_TH+1,NX_S_TH/2
            CU(I,K,J)=0.d0
          END DO
        END DO
      END DO 

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Z_TO_PHYSICAL_TH(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms along the z direction to Fourier space
C The input and output are real arrays with the input array packed
C to contain both the real and imagingary parts of the transformed array
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
      INTEGER  I, J, K
      COMPLEX*16   CU(0:NX_S_TH/2,0:NZ_TH+1,0:NY_S_TH)
      COMPLEX*16   U (0:NX_S_TH/2,0:NZ_TH+1,0:NY_S_TH)

      COMPLEX*16 CZX_PLANE(0:NZ_TH,0:NX_S_TH/2)

C First, unpack data into the CZX_PLANE complex temporary
C storage variable (dealiasing in the process)
C Then, perform a complex -> complex transform in the z-direction

      DO J=0,NY_S_TH
        DO I=0,NKX_S_TH
         DO K=0,NKZ_TH
            CZX_PLANE(K,I)=CU(I,K,J)
          END DO
          DO K=NKZ_TH+1,NZM_TH-NKZ_TH
            CZX_PLANE(K,I)=CMPLX(0.0,0.0)
          END DO
          DO K=1,NKZ_TH
            CZX_PLANE(NZM_TH-NKZ_TH+K,I)=CU(I,NKZ_TH+K,J)
          END DO
        END DO
        CALL FFTWND_F77(FFTW_Z_TO_P_PLAN_TH, NKX_S_TH+1,
     *    CZX_PLANE(0,0), 1, NZ_TH+1, CZX_PLANE(0,0), 1, NZ_TH+1)

        DO K=0,NZM_TH
          DO I=0,NKX_S_TH
            U(I,K,J)=CZX_PLANE(K,I)
          END DO
          DO I=NKX_S_TH+1,NX_S_TH/2
            U(I,K,J)=0.d0
          END DO
        END DO
       END DO

       RETURN
       END 

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Z_TO_FOURIER_TH(U,CU)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms in the z-direction
C The input and output should be in physical space with the output
C packed to hold the real and imaginary parts
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
      INTEGER I, J, K
      COMPLEX*16 U(0:NKX_S_TH,0:NZ_TH+1,0:NY_S_TH)
      COMPLEX*16 CU(0:NKX_S_TH,0:NZ_TH+1,0:NY_S_TH)

      COMPLEX*16 IN_PLANE(0:NZ_TH,0:NKX_S_TH)
      COMPLEX*16 OUT_PLANE(0:NZ_TH,0:NKX_S_TH)

C First, put the data into the CZX_PLANE temporary storage variable,
C perform a real -> complex transform in the z direction.


      DO J=0,NY_S_TH

        DO I=0,NKX_S_TH
          DO K=0,NZM_TH
            IN_PLANE(K,I)=U(I,K,J)
          END DO
        END DO


        CALL FFTWND_F77(FFTW_Z_TO_F_PLAN_TH, NKX_S_TH+1,
     *    IN_PLANE(0,0), 1, NZ_TH+1, IN_PLANE(0,0), 1, NZ_TH+1)


C Scale by NZ (necessary by FFTW convention)

        DO I=0,NKX_S_TH
          DO K=0,NKZ_TH
            CU(I,K,J)=IN_PLANE(K,I)/RNZ_TH
          END DO
          DO K=1,NKZ_TH
            CU(I,NKZ_TH+K,J)=IN_PLANE(NZM_TH-NKZ_TH+K,I)/RNZ_TH
          END DO
        END DO

      END DO

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Y_TO_FOURIER_TH(U,CU)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms in the y-direction
C The input and output should be in physical space with the output
C packed to hold the real and imaginary parts
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
      INTEGER I, J, K
      COMPLEX*16  U (0:NKX_S_TH,0:TNKZ_S_TH,0:NY_TH+1)
      COMPLEX*16  CU(0:NKX_S_TH,0:TNKZ_S_TH,0:NY_TH+1)

      COMPLEX*16 CYZ_PLANE(0:NY_TH,0:TNKZ_S_TH)

C First, put the data into the CYZ_PLANE temporary storage variable,
C perform a real -> complex transform in the y direction.

      DO I=0,NKX_S_TH
        DO K=0,TNKZ_S_TH
          DO J=0,NY_TH-1
            CYZ_PLANE(J,K)=U(I,K,J)
          END DO
        END DO

        CALL FFTWND_F77(FFTW_Y_TO_F_PLAN_TH, TNKZ_S_TH+1,
     *    CYZ_PLANE(0,0), 1, NY_TH+1, CYZ_PLANE(0,0), 1, NY_TH+1)

C Scale by NY (necessary by FFTW convention)
        DO K=0,TNKZ_S_TH
          DO J=0,NKY_TH
            CU(I,K,J)=CYZ_PLANE(J,K)/RNY_TH
          END DO
          DO J=1,NKY_TH
            CU(I,K,NKY_TH+J)=CYZ_PLANE(NYM_TH-NKY_TH+J,K)/RNY_TH
          END DO
        END DO            
      END DO

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Y_TO_PHYSICAL_TH(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms along the y direction to Fourier space
C The input and output are real arrays with the input array packed
C to contain both the real and imagingary parts of the transformed array
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
      INTEGER  I, J, K
      COMPLEX*16   CU(0:NX_S_TH/2,0:NZ_S_TH,0:NY_TH+1)
      COMPLEX*16   U (0:NX_S_TH/2,0:NZ_S_TH,0:NY_TH+1)

      COMPLEX*16 CYZ_PLANE(0:NY_TH,0:NZ_S_TH)

C First, unpack data into the CYZ_PLANE complex temporary
C storage variable (dealiasing in the process)
C Then, perform a complex -> real transform in the y-direction

      CYZ_PLANE(:,:)=(0.d0,0.d0)

      DO I=0,NKX_S_TH
        DO K=0,NZ_S_TH
          DO J=0,NKY_TH
            CYZ_PLANE(J,K)=CU(I,K,J)
          END DO
          DO J=NKY_TH+1,NYM_TH-NKY_TH
            CYZ_PLANE(J,K)=CMPLX(0.0,0.0)
          END DO
          DO J=1,NKY_TH
            CYZ_PLANE(NYM_TH-NKY_TH+J,K)=CU(I,K,NKY_TH+J)
          END DO
        END DO
        CALL FFTWND_F77(FFTW_Y_TO_P_PLAN_TH, NZ_S_TH+1,
     *    CYZ_PLANE(0,0), 1, NY_TH+1, CYZ_PLANE(0,0), 1, NY_TH+1)
        DO K=0,NZ_S_TH
          DO J=0,NYM_TH
            U(I,K,J)=CYZ_PLANE(J,K)
          END DO
        END DO
       END DO

       RETURN
       END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_X_TO_FOURIER_TH(U,CU)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
      INTEGER I, J, K
      REAL*8     U (0:NX_TH+1,0:NZ_S_TH,0:NY_S_TH)
      COMPLEX*16 CU(0:NX_TH/2,0:NZ_S_TH,0:NY_S_TH)

C Looping over the planes of interest, simply perform a real -> complex
C transform in place in the big storage array, scaling appropriately.

      DO J=0,NY_S_TH
      CALL RFFTWND_F77_REAL_TO_COMPLEX(FFTW_X_TO_F_PLAN_TH,(NZ_S_TH+1),
     *    U(0,0,J), 1, NX_TH+2, CU(0,0,J), 1, NX_TH/2+1)
        DO K=0,NZ_S_TH
          DO I=0,NKX_TH
            CU(I,K,J)=CU(I,K,J)/RNX_TH
          END DO
        END DO
      END DO

      DO j=0,NY_S_TH
        DO K=0,NZ_S_TH
          DO I=NKX_TH+1,NX_TH/2
            CU(I,K,J)=0.d0
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_X_TO_PHYSICAL_TH(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) to physical space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'

      INTEGER I, J, K
      REAL*8     U (0:NX_TH+1,0:NZ_S_TH,0:NY_S_TH)
      COMPLEX*16 CU(0:NX_TH/2,0:NZ_S_TH,0:NY_S_TH)

C Looping over the planes of interest, simply set the higher wavenumbers to
C zero and then perform a complex -> real transform in place in the big
C storage array.

      DO J=0,NY_S_TH
        DO K=0,NZ_S_TH
          DO I=NKX_TH+1,NX_TH/2
C Perform De-aliasing on the largest third of the wavenumbers
            CU(I,K,J)=0.
          END DO
        END DO

      CALL RFFTWND_F77_COMPLEX_TO_REAL(FFTW_X_TO_P_PLAN_TH,(NZ_S_TH+1),
     *    CU(0,0,J), 1, NX_TH/2+1, U(0,0,J), 1, NX_TH+2)

      END DO


      RETURN
      END


! The following subroutine is for taking a velocity field in Fourier space
! and interpolating it to the scalar (TH) grid in Physical space
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZY_MPI_TO_PHYSICAL_INTERP(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'

      INTEGER JMIN, JMAX, KMIN, KMAX, I, J, K
C Input array is in Fourier space and local in X on the velocity grid
      COMPLEX*16 CU(0:NX_S/2,0:NZ_S,0:NY+1)
C Output array is in physical space and local in Y on the TH grid
      REAL*8  U (0:NX_TH+1,0:NZ_S_TH,0:NY_S_TH+1)

C Intermediate arrays local X, and Z, These are all equivalenced
      COMPLEX*16 CU_X(0:NX/2,0:NZ_S_TH,0:NY_S_TH)
      COMPLEX*16 CU_X_TH(0:NX_TH/2,0:NZ_S_TH,0:NY_S_TH)
      REAL*8 U_X(0:NX_TH+1,0:NZ_S_TH,0:NY_S_TH)

C CU_Y and CU_Z can't be equivalenced since they are used
C in the transpose routines
      COMPLEX*16 CU_Y(0:NX_S/2,0:NZ_S,0:NY_TH+1)

      COMPLEX*16 CU_Z(0:NX_S/2,0:NZ+1,0:NY_S_TH)

      COMPLEX*16 CU_Z_TH(0:NX_S/2,0:NZ_TH+1,0:NY_S_TH)

C Equivalence the intermediate arrrays to avoid wasting memory
C The FFTs are done in-place, so this is safe 
      EQUIVALENCE(CU_X_TH,U_X)


C First, pad the velocity array with zeros
      CU_Y(:,:,:)=0.d0
      DO J=0,NKY
        CU_Y(:,:,J)=CU(:,:,J)
      END DO
      DO J=0,NKY-1
        CU_Y(:,:,TNKY_TH-J)=CU(:,:,TNKY-J)
      END DO
C Inverse transform in the x-direction:
      CALL FFT_Y_TO_PHYSICAL_INTERP(CU_Y,CU_Y)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C Perform a parallel transpose to get data stored locally in the z-direction
      CALL TRANSPOSE_MPI_Y_TO_Z_INTERP(CU_Y,CU_Z)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C Pad the high resolution array with zeros
      CU_Z_TH(:,:,:)=0.d0
      DO K=0,NKZ
        CU_Z_TH(:,K,:)=CU_Z(:,K,:)
      END DO
      DO K=0,NKZ-1
        CU_Z_TH(:,TNKZ_TH-K,:)=CU_Z(:,TNKZ-K,:)
      END DO

C Inverse transform in the z-direction:
      CALL FFT_Z_TO_PHYSICAL_INTERP(CU_Z_TH,CU_Z_TH)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C Perform a parallel transpose to get data stored locally in the y-direction
      CALL TRANSPOSE_MPI_Z_TO_X_INTERP(CU_Z_TH,CU_X)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C Pad the high resolution array with zeros
      CU_X_TH(:,:,:)=0.d0
      DO I=0,NKX
        CU_X_TH(I,:,:)=CU_X(I,:,:)
      END DO

C Inverse transform in the y-direction:
      CALL FFT_X_TO_PHYSICAL_TH(CU_X_TH,U_X)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C Now, transfer from U_X to the output array
      DO J=0,NY_S_TH
        DO K=0,NZ_S_TH
          DO I=0,NX_TH+1
            U(I,K,J)=U_X(I,K,J)
          END DO
        END DO
      END DO

C Set data above NXM=0
      DO I=NX_TH,NX_TH+1
        U(I,0:NZ_S,0:NY_S)=0.d0
      END DO
C Set data above NYM=0
      DO J=NY_TH-INT(RANK/NP_S)*(NY_S_TH+1),NY_S_TH+1
        U(0:NX_TH+1,0:NZ_S_TH,J)=0.d0
      END DO
C Set data above NZM=0
      DO K=NZ_TH-MOD(RANK,NP_S)*(NZ_S_TH+1),NZ_S_TH
        U(0:NX_TH+1,K,0:NY_S_TH)=0.d0
      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Z_TO_PHYSICAL_INTERP(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms along the z direction to Fourier space
C The input and output are real arrays with the input array packed
C to contain both the real and imagingary parts of the transformed array
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
      INTEGER  I, J, K
      COMPLEX*16   CU(0:NX_S/2,0:NZ_TH+1,0:NY_S_TH)
      COMPLEX*16   U (0:NX_S/2,0:NZ_TH+1,0:NY_S_TH)

      COMPLEX*16 CZX_PLANE(0:NZ_TH,0:NX_S/2)

C First, unpack data into the CZX_PLANE complex temporary
C storage variable (dealiasing in the process)
C Then, perform a complex -> complex transform in the z-direction

      DO J=0,NY_S_TH
        DO I=0,NKX_S
         DO K=0,NKZ_TH
            CZX_PLANE(K,I)=CU(I,K,J)
          END DO
          DO K=NKZ_TH+1,NZM_TH-NKZ_TH
            CZX_PLANE(K,I)=CMPLX(0.0,0.0)
          END DO
          DO K=1,NKZ_TH
            CZX_PLANE(NZM_TH-NKZ_TH+K,I)=CU(I,NKZ_TH+K,J)
          END DO
        END DO
        CALL FFTWND_F77(FFTW_Z_TO_P_PLAN_TH, NKX_S+1,
     *    CZX_PLANE(0,0), 1, NZ_TH+1, CZX_PLANE(0,0), 1, NZ_TH+1)
        DO K=0,NZM_TH
          DO I=0,NKX_S
            U(I,K,J)=CZX_PLANE(K,I)
          END DO
          DO I=NKX_S+1,NX_S/2
            U(I,K,J)=0.d0
          END DO
        END DO
       END DO

       RETURN
       END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_Y_TO_PHYSICAL_INTERP(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms along the y direction to Fourier space
C The input and output are real arrays with the input array packed
C to contain both the real and imagingary parts of the transformed array
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
      INTEGER  I, J, K
      COMPLEX*16   CU(0:NX_S/2,0:NZ_S,0:NY_TH+1)
      COMPLEX*16   U (0:NX_S/2,0:NZ_S,0:NY_TH+1)

      COMPLEX*16 CYZ_PLANE(0:NY_TH,0:NZ_S)

C First, unpack data into the CYZ_PLANE complex temporary
C storage variable (dealiasing in the process)
C Then, perform a complex -> real transform in the y-direction

      CYZ_PLANE(:,:)=(0.d0,0.d0)

      DO I=0,NKX_S
        DO K=0,NZ_S
          DO J=0,NKY_TH
            CYZ_PLANE(J,K)=CU(I,K,J)
          END DO
          DO J=NKY_TH+1,NYM_TH-NKY_TH
            CYZ_PLANE(J,K)=CMPLX(0.0,0.0)
          END DO
          DO J=1,NKY_TH
            CYZ_PLANE(NYM_TH-NKY_TH+J,K)=CU(I,K,NKY_TH+J)
          END DO
        END DO
        CALL FFTWND_F77(FFTW_Y_TO_P_PLAN_TH, NZ_S+1,
     *    CYZ_PLANE(0,0), 1, NY_TH+1, CYZ_PLANE(0,0), 1, NY_TH+1)
        DO K=0,NZ_S
          DO J=0,NYM_TH
            U(I,K,J)=CYZ_PLANE(J,K)
          END DO
        END DO
       END DO

       RETURN
       END





C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZ_TO_FOURIER_TH(U,CU,JMIN,JMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 2 directions) planes JMIN-JMAX to Fourier space.

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZ_TO_PHYSICAL_TH(CU,U,JMIN,JMAX)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 2 directions) planes JMIN-JMAX to physical space.
      RETURN
      END
C******************************************************************************|
C--------> The transform routines for the fully-periodic box follow. <---------|
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZY_TO_FOURIER_TH(U,CU)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 3 directions) the entire box to Fourier space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FFT_XZY_TO_PHYSICAL_TH(CU,U)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This routine transforms (in 3 directions) the entire box to physical space.
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'


      RETURN
      END





