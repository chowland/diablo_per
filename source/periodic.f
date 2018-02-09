
C periodic.f, the fully-periodic-box solvers for diablo.           VERSION 0.9
C
C These solvers were written by Tom Bewley and John Taylor.
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE 'header'

!      IF (N_TH.gt.0) THEN
!      EPSILON_TARGET=((1.d0/(DX(1)*0.25d0))**4.d0)
!     &              *(NU**3.d0)*(700.d0)**(-2.d0)
!      EPSILON_TARGET=((1.d0/(2.d0*DX_TH(1)))**4.d0)
!     &              *(NU**3.d0)*(PR(1)/1.d0)**(-2.d0)
! Re_lambda=35:
      ETA_TARGET=LX/92.d0
! Re_lambda=25:
      ETA_TARGET=0.111d0
! High Reynolds number, SC=30
      ETA_TARGET=0.0336078d0
      if (rank.eq.0) write(*,*) 'TARGET KOLMOGOROV SCALE: ',ETA_TARGET
!      EPSILON_TARGET=((3.d0/(2.d0*DX(1)))**4.d0)
!     &        *(NU**3.d0)*(100)**(-2.d0)
!      write(*,*) 'EPSILON_TARGET: ',EPSILON_TARGET
!      END IF

C Define variables for the Geophysical case
      PI=4.D0*ATAN(1.D0)
      PHI=2.D0*PI*PHI/360.D0
      GAMMA=2.D0*PI*GAMMA/360.D0
      C_SIN=COS(PHI)*SIN(GAMMA)/SIN(PHI)
      C_COS=COS(PHI)*COS(GAMMA)/SIN(PHI)

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_PER_1
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Main time-stepping algorithm for the fully periodic case
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE 'header'

      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6
      INTEGER I,J,K,N

C Start with data local in the X-direction

C Compute the RHS in Fourier Space, CRi.

C First, define the constants used for time-stepping
      TEMP1=NU*H_BAR(RK_STEP)/2.0
      TEMP2=BETA_BAR(RK_STEP)*H_BAR(RK_STEP)
      TEMP3=ZETA_BAR(RK_STEP)*H_BAR(RK_STEP)
      TEMP4=H_BAR(RK_STEP)

      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
C Start with the explicit part of the Crank-Nicolson viscous term and
C  the pressure gradient treated with Explicit Euler:
            TEMP5=1-TEMP1*(KX2_S(I)+KY2(J)+KZ2_S(K))**BETA

            CR1(I,K,J)=TEMP5*CU1(I,K,J)-TEMP4*(CIKX_S(I)*CP(I,K,J))
            CR2(I,K,J)=TEMP5*CU2(I,K,J)-TEMP4*(CIKY(J)*CP(I,K,J))
            CR3(I,K,J)=TEMP5*CU3(I,K,J)-TEMP4*(CIKZ_S(K)*CP(I,K,J))
          END DO
        END DO
      END DO

C For each scalar, start with the explict part of the Crank-Nicolson
C diffusive term for each scalar
      DO N=1,N_TH
        DO J=0,TNKY_TH
          DO K=0,TNKZ_S_TH
            DO I=0,NKX_S_TH
            TEMP6=1-(TEMP1/PR(N))*(KX2_S_TH(I)
     &               +KY2_TH(J)+KZ2_S_TH(K))**BETA
     &               -REACTION(N)*TEMP1
              CRTH(I,K,J,N)=TEMP6*CTH(I,K,J,N)
            END DO
          END DO
        END DO
      END DO

      IF (RK_STEP .GT. 1) THEN
        DO J=0,TNKY
          DO K=0,TNKZ_S
            DO I=0,NKX_S
C Add the term: ZETA_BAR(RK_STEP)*R(U(RK_STEP-1))
              CR1(I,K,J)=CR1(I,K,J)+TEMP3*CF1(I,K,J)
              CR2(I,K,J)=CR2(I,K,J)+TEMP3*CF2(I,K,J)
              CR3(I,K,J)=CR3(I,K,J)+TEMP3*CF3(I,K,J)
            END DO
          END DO
        END DO
        DO N=1,N_TH
        DO J=0,TNKY_TH
          DO K=0,TNKZ_S_TH
            DO I=0,NKX_S_TH
                CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP3*CFTH(I,K,J,N)
              END DO
            END DO
          END DO
        END DO
      END IF


C If we are considering a linear background scalar gradient then add
C the term owing to advection of the background state.
C This allows us to consider a mean scalar gradient (ie stratification)
C even though the vertical boundary conditions are periodic.
C (In this case the passive scalar is a perturbation from a linear
C gradient. This gradient and the vertical domain size are used to
C make the passive scalar nondimensional, so here the nondimensional
C gradient is equal to one
      DO N=1,N_TH
      IF (BACKGROUND_GRAD(N)) THEN
C If there is a background scalar gradient add advection term:
C Start by setting CFTH to zero:
      CFTH(:,:,:,N)=0.d0
C Now, add the velocity advection only to velocity modes
      DO J=0,TNKY_TH
        DO K=0,TNKZ_S_TH
          DO I=0,NKX_S_TH
            CFTH(I,K,J,N)=-CU2(I,K,J)
          END DO
        END DO
      END DO
      ELSE
! Otherwise don't
        CFTH(:,:,:,N)=0.d0
      END IF
      END DO


C For Geophysical applications, add the Coriolis term
C Note, on an f-plane under the traditional approximation, C_SIN=C_COS=0
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
           CF1(I,K,J)=I_RO_TAU*(C_SIN*CU2(I,K,J)+CU3(I,K,J))
           CF2(I,K,J)=I_RO_TAU*(-1.d0*C_COS*CU3(I,K,J)+C_SIN*CU1(I,K,J))
           CF3(I,K,J)=I_RO_TAU*(C_COS*CU2(I,K,J)-CU1(I,K,J))
!           CF1(I,K,J)=I_RO_TAU*CU3(I,K,J)
!           CF2(I,K,J)=0.d0
!           CF3(I,K,J)=-1.d0*I_RO_TAU*CU1(I,K,J)
          END DO
        END DO
      END DO

!      call USER_RHS_PER_FOURIER


C Transform the scalar concentration to physical space
      DO N=1,N_TH
        IF (USE_MPI) THEN
          CALL FFT_XZY_MPI_TO_PHYSICAL_TH(CTH(0,0,0,N),TH(0,0,0,N))
        ELSE
          CALL FFT_XZY_TO_PHYSICAL_TH(CTH(0,0,0,N),TH(0,0,0,N))
        END IF
      END DO


C Compute the nonlinear terms for the passive scalar equation
C Do this before the nonlinear momentum terms to use Fi as a working
C  array before using it for the momentum equation.

      DO N=1,N_TH
      CALL FFT_XZY_MPI_TO_PHYSICAL_INTERP(CU1,STH1)
        DO J=0,NY_S_TH
          DO K=0,NZ_S_TH
            DO I=0,NXM_TH
              STH1(I,K,J)=STH1(I,K,J)*TH(I,K,J,N)
            END DO
          END DO
        END DO
        IF (USE_MPI) THEN
          CALL FFT_XZY_MPI_TO_FOURIER_TH(STH1,CSTH1)
        ELSE
          CALL FFT_XZY_TO_FOURIER_TH(STH1,CSTH1)
        END IF
        DO J=0,TNKY_TH
          DO K=0,TNKZ_S_TH
            DO I=0,NKX_S_TH
              CFTH(I,K,J,N)=CFTH(I,K,J,N)-CIKX_S_TH(I)*CSTH1(I,K,J)
            END DO
          END DO
        END DO

        CALL FFT_XZY_MPI_TO_PHYSICAL_INTERP(CU2,STH1)
        DO J=0,NY_S_TH
          DO K=0,NZ_S_TH
            DO I=0,NXM_TH
              STH1(I,K,J)=STH1(I,K,J)*TH(I,K,J,N)
            END DO
          END DO
       END DO
        IF (USE_MPI) THEN
          CALL FFT_XZY_MPI_TO_FOURIER_TH(STH1,CSTH1)
        ELSE
          CALL FFT_XZY_TO_FOURIER_TH(STH1,CSTH1)
        END IF
        DO J=0,TNKY_TH
          DO K=0,TNKZ_S_TH
            DO I=0,NKX_S_TH
              CFTH(I,K,J,N)=CFTH(I,K,J,N)-CIKY_TH(J)*CSTH1(I,K,J)
            END DO
          END DO
        END DO

        CALL FFT_XZY_MPI_TO_PHYSICAL_INTERP(CU3,STH1)
        DO J=0,NY_S_TH
          DO K=0,NZ_S_TH
            DO I=0,NXM_TH
              STH1(I,K,J)=STH1(I,K,J)*TH(I,K,J,N)
            END DO
          END DO
        END DO
        IF (USE_MPI) THEN
          CALL FFT_XZY_MPI_TO_FOURIER_TH(STH1,CSTH1)
        ELSE
          CALL FFT_XZY_TO_FOURIER_TH(STH1,CSTH1)
        END IF
        DO J=0,TNKY_TH
          DO K=0,TNKZ_S_TH
            DO I=0,NKX_S_TH
              CFTH(I,K,J,N)=CFTH(I,K,J,N)-CIKZ_S_TH(K)*CSTH1(I,K,J)
            END DO
          END DO
        END DO
      END DO

C Add R-K terms for the TH equation to the RHS
      DO N=1,N_TH
        DO J=0,TNKY_TH
          DO K=0,TNKZ_S_TH
            DO I=0,NKX_S_TH
              CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP2*CFTH(I,K,J,N)
            END DO
          END DO
        END DO
      END DO
C The RHS vector for the TH equation is now ready

C Inverse transform to physical space to compute the nonlinear terms
      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_PHYSICAL(CU1,U1)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CU2,U2)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CU3,U3)
      ELSE
        CALL FFT_XZY_TO_PHYSICAL(CU1,U1)
        CALL FFT_XZY_TO_PHYSICAL(CU2,U2)
        CALL FFT_XZY_TO_PHYSICAL(CU3,U3)
      END IF

C Compute the nonlinear terms for the momentum equations
        DO J=0,NY_S
          DO K=0,NZ_S
            DO I=0,NXM
              S1(I,K,J)=U1(I,K,J)*U1(I,K,J)
            END DO
          END DO
        END DO
        IF (USE_MPI) THEN
          CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
        ELSE
          CALL FFT_XZY_TO_FOURIER(S1,CS1)
        END IF
        DO J=0,TNKY
          DO K=0,TNKZ_S
            DO I=0,NKX_S
              CF1(I,K,J)=CF1(I,K,J)-CIKX_S(I)*CS1(I,K,J)
            END DO
          END DO
        END DO

        DO J=0,NY_S
          DO K=0,NZ_S
            DO I=0,NXM
              S1(I,K,J)=U1(I,K,J)*U2(I,K,J)
            END DO
          END DO
        END DO
        IF (USE_MPI) THEN
          CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
        ELSE
          CALL FFT_XZY_TO_FOURIER(S1,CS1)
        END IF
        DO J=0,TNKY
          DO K=0,TNKZ_S
            DO I=0,NKX_S
              CF1(I,K,J)=CF1(I,K,J)-CIKY(J)*CS1(I,K,J)
              CF2(I,K,J)=CF2(I,K,J)-CIKX_S(I)*CS1(I,K,J)
            END DO
          END DO
        END DO

        DO J=0,NY_S
          DO K=0,NZ_S
            DO I=0,NXM
              S1(I,K,J)=U1(I,K,J)*U3(I,K,J)
            END DO
          END DO
        END DO
        IF (USE_MPI) THEN
          CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
        ELSE
          CALL FFT_XZY_TO_FOURIER(S1,CS1)
        END IF
        DO J=0,TNKY
          DO K=0,TNKZ_S
            DO I=0,NKX_S
              CF1(I,K,J)=CF1(I,K,J)-CIKZ_S(K)*CS1(I,K,J)
              CF3(I,K,J)=CF3(I,K,J)-CIKX_S(I)*CS1(I,K,J)
            END DO
          END DO
        END DO

        DO J=0,NY_S
          DO K=0,NZ_S
            DO I=0,NXM
              S1(I,K,J)=U2(I,K,J)*U2(I,K,J)
            END DO
          END DO
        END DO
        IF (USE_MPI) THEN
          CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
        ELSE
          CALL FFT_XZY_TO_FOURIER(S1,CS1)
        END IF
        DO J=0,TNKY
          DO K=0,TNKZ_S
            DO I=0,NKX_S
              CF2(I,K,J)=CF2(I,K,J)-CIKY(J)*CS1(I,K,J)
            END DO
          END DO
        END DO

       DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            S1(I,K,J)=U2(I,K,J)*U3(I,K,J)
          END DO
        END DO
       END DO
       IF (USE_MPI) THEN
         CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
       ELSE
         CALL FFT_XZY_TO_FOURIER(S1,CS1)
       END IF
       DO J=0,TNKY
         DO K=0,TNKZ_S
          DO I=0,NKX_S
            CF2(I,K,J)=CF2(I,K,J)-CIKZ_S(K)*CS1(I,K,J)
            CF3(I,K,J)=CF3(I,K,J)-CIKY(J)*CS1(I,K,J)
          END DO
         END DO
       END DO

       DO J=0,NY_S
         DO K=0,NZ_S
          DO I=0,NXM
            S1(I,K,J)=U3(I,K,J)*U3(I,K,J)
          END DO
         END DO
       END DO
       IF (USE_MPI) THEN
         CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
       ELSE
         CALL FFT_XZY_TO_FOURIER(S1,CS1)
       END IF
       DO J=0,TNKY
         DO K=0,TNKZ_S
          DO I=0,NKX_S
            CF3(I,K,J)=CF3(I,K,J)-CIKZ_S(K)*CS1(I,K,J)
          END DO
         END DO
       END DO
C Done with the computation of nonlinear terms

C If the scalar is active (RI_TAU NE 0), add the bouyancy forcing term
C  as explicit R-K
       DO N=1,N_TH
       IF (RI_TAU(N).NE.0.d0) THEN
C Fisrt, convert back to Fourier space
       IF (USE_MPI) THEN
         CALL FFT_XZY_MPI_TO_FOURIER_TH(TH(0,0,0,N),CTH(0,0,0,N))
       ELSE
         CALL FFT_XZY_TO_FOURIER_TH(TH(0,0,0,N),CTH(0,0,0,N))
       END IF
C Add only to the velocity modes:
C ASSUME that the fine scale scalar structure does not affect the velocity
! Note, this needs to be modified before it can be used with varying grid sizes
       DO J=0,TNKY
         DO K=0,TNKZ_S
           DO I=0,NKX_S
             CF2(I,K,J)=CF2(I,K,J)+RI_TAU(N)*CTH(I,K,J,N)
           END DO
         END DO
       END DO
       END IF
       END DO

       !call USER_RHS_PER_FOURIER

C Add some forcing to the system to keep the Batchelor scale fixed
!      EK=0.d0
!      DO J=0,NY_S
!        DO K=0,NZ_S
!          DO I=0,NXM
!              EK=EK+U1(I,K,J)**2.d0+U2(I,K,J)**2.d0+U3(I,K,J)**2.d0
!          END DO
!        END DO
!      END DO
C Note, that each cell has the same volume, so we can just average over all points
!      EK=EK/dble(NX*NY*NZ)
!      IF (USE_MPI) THEN
!        CALL MPI_ALLREDUCE(EK,EK_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM
!     &                     ,MPI_COMM_WORLD,IERROR)
!      END IF
!      RE_LAMBDA=EK*SQRT(15.d0/(EPSILON_TARGET*NU))
!      write(*,*) 'RE_LAMBDA: ',RE_LAMBDA
!      IF (RE_LAMBDA.lt.35.d0) THEN
!      EPSILON_TARGET=EPSILON_TARGET+0.0002d0*ABS(35.d0-RE_LAMBDA)/35.d0
!      ELSE IF (RE_LAMBDA.gt.35.d0) THEN
!      EPSILON_TARGET=EPSILON_TARGET-0.0002d0*ABS(35.d0-RE_LAMBDA)/35.d0
!      END IF
!      WRITE(*,*) 'EPSILON_TARGET: ',EPSILON_TARGET
!      EPSILON_TARGET=15.d0*EK**2.d0/(35.d0**2.d0*NU)

! Scale EK by an amount to compensate for dissipation from 2/3 de-aliasing:
!      EK=0.8d0*EK_sum

!      DO J=0,NY_S
!        DO K=0,NZ_S
!          DO I=0,NXM
!            S1(I,K,J)=(EPSILON_TARGET/EK)*U1(I,K,J)
!          END DO
!        END DO
!      END DO
!      IF (USE_MPI) THEN
!        CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
!      ELSE
!        CALL FFT_XZY_TO_FOURIER(S1,CS1)
!      END IF
!      DO J=0,TNKY
!        DO K=0,TNKZ_S
!          DO I=0,NKX_S
!            CF1(I,K,J)=CF1(I,K,J)+CS1(I,K,J)
!          END DO
!        END DO
!      END DO
!      DO J=0,NY_S
!       DO K=0,NZ_S
!          DO I=0,NXM
!            S1(I,K,J)=(EPSILON_TARGET/EK)*U2(I,K,J)
!          END DO
!        END DO
!      END DO
!      IF (USE_MPI) THEN
!        CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
!      ELSE
!        CALL FFT_XZY_TO_FOURIER(S1,CS1)
!      END IF
!      DO J=0,TNKY
!        DO K=0,TNKZ_S
!          DO I=0,NKX_S
!            CF2(I,K,J)=CF2(I,K,J)+CS1(I,K,J)
!          END DO
!        END DO
!      END DO
!      DO J=0,NY_S
!        DO K=0,NZ_S
!          DO I=0,NXM
!            S1(I,K,J)=(EPSILON_TARGET/EK)*U3(I,K,J)
!          END DO
!        END DO
!      END DO
!      IF (USE_MPI) THEN
!        CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
!      ELSE
!        CALL FFT_XZY_TO_FOURIER(S1,CS1)
!      END IF
!      DO J=0,TNKY
!        DO K=0,TNKZ_S
!          DO I=0,NKX_S
!            CF3(I,K,J)=CF3(I,K,J)+CS1(I,K,J)
!          END DO
!        END DO
!      END DO

C Now, add the R-K terms to the RHS
C Note, this has already been done for the scalar field TH
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CR1(I,K,J)=CR1(I,K,J)+TEMP2*CF1(I,K,J)
            CR2(I,K,J)=CR2(I,K,J)+TEMP2*CF2(I,K,J)
            CR3(I,K,J)=CR3(I,K,J)+TEMP2*CF3(I,K,J)
          END DO
        END DO
      END DO
C Computation of CRi complete.

C If Variable timestepping and done with one full R-K step, update
C DELTA_T based on the specified CFL number
! Note, this change will not take effect until the next timestep
! since the TEMP variables have already been defined
      IF ((VARIABLE_DT).and.(RK_STEP.eq.3)
     &        .and.(MOD(TIME_STEP,UPDATE_DT).EQ.0)) THEN
        IF (.NOT.USE_MPI) THEN
          CALL COURANT
        ELSE
          CALL COURANT_MPI
        END IF
      END IF

      CU1=(0.d0,0.d0)
      CU2=(0.d0,0.d0)
      CU3=(0.d0,0.d0)
      CTH=(0.d0,0.d0)

C Now solve the implicit system for the intermediate field.
C (In the fully-periodic case, this is easy!)
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            TEMP5=1.d0+TEMP1*(KX2_S(I)+KY2(J)+KZ2_S(K))**BETA
            CU1(I,K,J)=CR1(I,K,J)/TEMP5
            CU2(I,K,J)=CR2(I,K,J)/TEMP5
            CU3(I,K,J)=CR3(I,K,J)/TEMP5
          END DO
        END DO
      END DO
      DO N=1,N_TH
        DO J=0,TNKY_TH
          DO K=0,TNKZ_S_TH
            DO I=0,NKX_S_TH
            TEMP6=1+(TEMP1/PR(N))*(KX2_S_TH(I)
     &             +KY2_TH(J)+KZ2_S_TH(K))**BETA
     &            + REACTION(N)*TEMP1
             CTH(I,K,J,N)=CRTH(I,K,J,N)/TEMP6
            END DO
          END DO
        END DO
      END DO
C First step of the Fractional Step algorithm complete.


C Begin second step of the Fractional Step algorithm, making u divergence free.

C Compute varphi, store in the variable CR1, and project velocity field
      CALL REM_DIV_PER

C Then, update the pressure with phi.
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CP(I,K,J)=CP(I,K,J)+CR1(I,K,J)/TEMP4
          END DO
        END DO
      END DO
C Second step of the Fractional Step algorithm complete.


      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_PER_2
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Alternative time-stepping algorithm for the fully periodic case.
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K
      REAL*8  TEMP5

C Compute phi, store in the variable CR1.
C Note the coefficient H_BAR is absorbed into phi.

      CR1=(0.d0,0.d0)

      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            TEMP5=-(KX2_S(I)+KY2(J)+KZ2_S(K)+EPS)
            CR1(I,K,J)=(CIKX_S(I)*CU1(I,K,J)+CIKY(J)*CU2(I,K,J)+
     *                  CIKZ_S(K)*CU3(I,K,J))/TEMP5
          END DO
        END DO
C Then update CUi to make velocity field divergence-free.
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CU1(I,K,J)=CU1(I,K,J)-CIKX_S(I)*CR1(I,K,J)
            CU2(I,K,J)=CU2(I,K,J)-CIKY(J)*CR1(I,K,J)
            CU3(I,K,J)=CU3(I,K,J)-CIKZ_S(K)*CR1(I,K,J)
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE POISSON_P_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K


      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CR1(I,K,J)=CIKX_S(I)*CU1(I,K,J)
            CR2(I,K,J)=CIKZ_S(K)*CU3(I,K,J)
          END DO
        END DO
      END DO

      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_PHYSICAL(CR1,R1)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CR2,R2)
      ELSE
        CALL FFT_XZY_TO_PHYSICAL(CR1,R1)
        CALL FFT_XZY_TO_PHYSICAL(CR2,R2)
      END IF

      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            P(I,K,J)=R1(I,K,J)*R1(I,K,J)+R1(I,K,J)*R2(I,K,J)+
     *               R2(I,K,J)*R2(I,K,J)
          END DO
        END DO
      END DO
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CF1(I,K,J)=CIKY(J)*CU1(I,K,J)
            CF2(I,K,J)=CIKX_S(I)*CU2(I,K,J)
            CF3(I,K,J)=CIKZ_S(K)*CU1(I,K,J)
            CR1(I,K,J)=CIKX_S(I)*CU3(I,K,J)
            CR2(I,K,J)=CIKZ_S(K)*CU2(I,K,J)
            CR3(I,K,J)=CIKY(J)*CU3(I,K,J)
          END DO
        END DO
      END DO
      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_PHYSICAL(CF1,F1)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CF2,F2)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CF3,F3)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CR1,R1)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CR2,R2)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CR3,R3)
      ELSE
        CALL FFT_XZY_TO_PHYSICAL(CF1,F1)
        CALL FFT_XZY_TO_PHYSICAL(CF2,F2)
        CALL FFT_XZY_TO_PHYSICAL(CF3,F3)
        CALL FFT_XZY_TO_PHYSICAL(CR1,R1)
        CALL FFT_XZY_TO_PHYSICAL(CR2,R2)
        CALL FFT_XZY_TO_PHYSICAL(CR3,R3)
      END IF
      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            P(I,K,J)=2*(P(I,K,J)+F1(I,K,J)*F2(I,K,J)
     *                          +F3(I,K,J)*R1(I,K,J)
     *                          +R2(I,K,J)*R3(I,K,J))
          END DO
        END DO
      END DO
      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_FOURIER(P,CP)
      ELSE
        CALL FFT_XZY_TO_FOURIER(P,CP)
      END IF

      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CP(I,K,J)=CP(I,K,J)/(KX2_S(I)+KY2(J)+KZ2_S(K)+EPS)
          END DO
        END DO
      END DO

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I, J, K
      REAL*8 RNUM1,RNUM2,RNUM3
      REAL*8 K0,K_MAG
      CHARACTER*60 FNAME
      REAL*8 b3,f03,E3,Sigma3,A3,alpha


C For an initial vortex, define the location of the centerline
      REAL*8 XC(0:NY+1),ZC(0:NY+1)

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

      IF (FLAVOR .EQ. 'Basic') THEN

      WRITE(6,*) 'Creating new flow from scratch.'

C Initialize random number generator
      CALL RANDOM_SEED

      IF (IC_TYPE.eq.0) THEN
C Initizlize the flow using a Taylor-Green vortex
C Nondimensionalize with U0 and 1/kappa

      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
!            IF (((J+INT(RANK/NP_S)*(NY_S+1)).le.NY)
!     &          .and.((K+MOD(RANK,NP_S)*(NZ_S+1)).le.NZ)) then
! Add a random phase
!            U1(I,K,J)=cos(2*pi*(GY_S(J))/LY)
!     &               *cos(2*pi*(GX(I))/LX)
!     &               *SIN(2*pi*(GZ_S(K))/LZ)
!            U2(I,K,J)=0.d0
!            U3(I,K,J)=-cos(2*pi*(GY_S(J))/LY)
!     &               *sin(2*pi*(GX(I))/LX)
!     &               *COS(2*pi*(GZ_S(K))/LZ)
!            else
             U1(i,k,j)=0.d0
             U2(i,k,j)=0.d0
             U3(i,k,j)=0.d0
!            end if
          END DO
        END DO
      END DO
      ELSE IF (IC_TYPE.eq.1) THEN
C Start with an ideal vortex centered in the domain
      DO J=0,NY_S
!        XC(J)=LX/2.+(LX/5.)*sin(2*PI*GY(J)/LY)
!        ZC(J)=LZ/2.+(LZ/5.)*sin(2*PI*GY(J)/LY)
        XC(J)=LX/2.
        ZC(J)=LZ/2.
        DO K=0,NZ_S
          DO I=0,NXM
            IF ((GX(I)-XC(j))**2.+(GZ_S(K)-ZC(j))**2..gt.0.1) then
! If we aren't too close to the vortex center
              U1(I,K,J)=-1.d0*(GZ_S(K)-ZC(j))
     &                /((GX(I)-XC(j))**2.+(GZ_S(K)-ZC(j))**2.)
              U3(I,K,J)=1.d0*(GX(I)-XC(j))
     &                /((GX(I)-XC(j))**2.+(GZ_S(K)-ZC(j))**2.)
              U2(I,K,J)=0.d0
            ELSE
! Otherwise:
              U1(I,K,J)=-1.d0*(GZ_S(K)-ZC(j))
     &                /0.1
              U3(I,K,J)=1.d0*(GX(I)-XC(j))
     &                /0.1
              U2(I,K,J)=0.d0
            END IF
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            CALL RANDOM_NUMBER(RNUM3)
            U1(I,K,J)=U1(I,K,J)+(RNUM1-0.5)*KICK
            U2(I,K,J)=U2(I,K,J)+(RNUM1-0.5)*KICK
            U3(I,K,J)=U3(I,K,J)+(RNUM1-0.5)*KICK
          END DO
        END DO
      END DO
      ELSE IF (IC_TYPE.eq.2) THEN
! Initilize with a plane wave
        theta0=45.d0*2.d0*PI/360.d0
        a0=0.85d0
        DO J=0,NY_S
          DO K=0,NZ_S
            DO I=0,NX+1
              phi0=GX(I)+GY_S(J)*tan(theta0)
              U1(I,K,J)=-a0*sin(theta0)*cos(phi0)
              U2(I,K,J)=a0*cos(theta0)*cos(phi0)
              U3(I,K,J)=0.d0
            END DO
          END DO
        END DO
      ELSE IF (IC_TYPE.eq.3) THEN
        ! Initialize with a GM spectrum of internal waves
        b3=65.
        f03=0.014*sqrt(RI_TAU(1))
        E3=6.3e-5
        Sigma3=0.468
        A3=sqrt(RI_TAU(1)*b3*E3*f03*sqrt(RI_TAU(1)-f03**2)*LX*LY*LZ)
     &          /sqrt(2*PI*Sigma3)
        do j=0,TNKY
          do k=0,TNKZ_S
            do i=0,NKX_S
              if ( (KX2_S(i)+KZ2_S(k)+KY2(j).le.100.) .and.
     &              (KY(j).ne.0) .and. (KX2_S(i)+KZ2_S(k).ne.0) ) then
                call RANDOM_NUMBER(alpha)
                alpha=2.*pi*alpha ! Random phase of each wave
                CS1(i,k,j)=A3*cexp(cmplx(0,alpha))*KY(j)/
     &         (sqrt(KX2_S(i)+KZ2_S(k))*
     &          (KX2_S(i)+KZ2_S(k)+KY2(j))**(1./4.)
     &         *(RI_TAU(1)*(KX2_S(i)+KZ2_S(k))+f03**2*KY2(j))
     &         *(KY2(j)+pi**2*9/b3**2))**(0.5)
                CU1(i,k,j)=CS1(i,k,j)*KY(j)*KX_S(i)/sqrt(
     &            (KX2_S(i)+KZ2_S(k)+KY2(j))*(KX2_S(i)+KZ2_S(k)))
                CU2(i,k,j)=CS1(i,k,j)*sqrt(KX2_S(i)+KZ2_S(k))/sqrt(
     &            KX2_S(i)+KZ2_S(k)+KY2(j))
                CU3(i,k,j)=CS1(i,k,j)*KY(j)*KZ_S(k)/sqrt(
     &            (KX2_S(i)+KZ2_S(k)+KY2(j))*(KX2_S(i)+KZ2_S(k)))
                CTH(i,k,j,1)=CS1(i,k,j)*CI/sqrt(RI_TAU(1))
              end if
            end do
          end do
        end do
      ELSE
        WRITE(*,*) 'Warning, Undefined Initial conditions in periodic.f'
      END IF

      if (IC_TYPE.ne.3) then
        IF (USE_MPI) THEN
          CALL FFT_XZY_MPI_TO_FOURIER(U1,CU1)
          CALL FFT_XZY_MPI_TO_FOURIER(U3,CU3)
          CALL FFT_XZY_MPI_TO_FOURIER(U2,CU2)
        ELSE
          CALL FFT_XZY_TO_FOURIER(U1,CU1)
          CALL FFT_XZY_TO_FOURIER(U2,CU2)
          CALL FFT_XZY_TO_FOURIER(U3,CU3)
        END IF
      end if

      IF (RANK.eq.0) THEN
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            RNUM1=(RNUM1-0.5d0)*2.d0
            RNUM2=(RNUM2-0.5d0)*2.d0
            CU1(I,K,J)=CU1(I,K,J)+CMPLX(RNUM1,RNUM2)*KICK
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            RNUM1=(RNUM1-0.5d0)*2.d0
            RNUM2=(RNUM2-0.5d0)*2.d0
            CU2(I,K,J)=CU2(I,K,J)+CMPLX(RNUM1,RNUM2)*KICK
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            RNUM1=(RNUM1-0.5d0)*2.d0
            RNUM2=(RNUM2-0.5d0)*2.d0
            CU3(I,K,J)=CU3(I,K,J)+CMPLX(RNUM1,RNUM2)*KICK
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            RNUM1=(RNUM1-0.5d0)*2.d0
            RNUM2=(RNUM2-0.5d0)*2.d0
          end do
        end do
      end do
      END IF

       ELSE IF (FLAVOR .EQ. 'Ensemble') THEN
       CALL CREATE_FLOW_ENSEM
       ELSE
       write(6,*) 'Unknown flavor, flow-field not created'
       endif

!      CALL SAVE_STATS_PER(.FALSE.)

        IF (USE_MPI) THEN
           CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        END IF

      if (RANK.eq.0) write(*,*) 'calling rem_div...'
      CALL REM_DIV_PER

!      CALL POISSON_P_PER

!      CALL SAVE_STATS_PER(.FALSE.)

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K,N

C Note, Since stratification is not permitted in the periodic flow field
C Any background stratification must be added to the governing equations

      IF ((IC_TYPE.eq.0).or.(IC_TYPE.eq.1).or.(IC_TYPE.eq.3)) then
        DO N=1,N_TH
          IF (CREATE_NEW_TH(N)) THEN
            DO J=0,NY_S_TH
              DO K=0,NZ_S_TH
                DO I=0,NXM_TH
                  TH(I,K,J,n)=0.d0
                END DO
              END DO
            END DO
          END IF
        END DO
      ELSE IF (IC_TYPE.eq.2) then
! for plane wave
        DO N=1,N_TH
          IF (CREATE_NEW_TH(N)) THEN
            theta0=45.d0*2.d0*PI/360.d0
            a0=0.85d0
            DO J=0,NY_S
              DO K=0,NZ_S
                DO I=0,NX+1
                  phi0=GX(I)+GY_S(J)*tan(theta0)
                  U1(I,K,J)=-a0*sin(theta0)*cos(phi0)
                  U2(I,K,J)=a0*cos(theta0)*cos(phi0)
                  U3(I,K,J)=0.d0
                  TH(I,K,J,1)=-a0*sin(phi0)
                END DO
              END DO
            END DO
          END IF
        END DO
       ELSE
         WRITE(*,*) 'UNKNOWN IC_TYPE IN CREATE_TH_PER'
       END IF

! Transfer to Fourier space
       DO N=1,N_TH
       IF (USE_MPI) THEN
         CALL FFT_XZY_MPI_TO_FOURIER_TH(TH(0,0,0,N),CTH(0,0,0,N))
       ELSE
         CALL FFT_XZY_TO_FOURIER_TH(TH(0,0,0,N),CTH(0,0,0,N))
       END IF
       END DO

       RETURN
       END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_GRID_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K

!      IF ((FLAVOR.NE.'Ensemble') .AND. (RANK.eq.0)) WRITE (6,*)
!     &            'Fourier in X'
      DO I=0,NX
        GX(I)=(I*LX)/NX
        DX(I)=LX/NX
        IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
      END DO
!      IF ((FLAVOR.NE.'Ensemble') .AND. (RANK.eq.0)) WRITE (6,*)
!     &            'Fourier in Z'
      DO K=0,NZ
        GZ(K)=(K*LZ)/NZ
        DZ(K)=LZ/NZ
        IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ(',K,') = ',GZ(K)
      END DO
!      IF ((FLAVOR.NE.'Ensemble') .AND. (RANK.eq.0)) WRITE (6,*)
!     &            'Fourier in Y'
      DO J=0,NY
        GY(J)=(J*LY)/NY
        DY(J)=LY/NY
        IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY(',J,') = ',GY(J)
      END DO

         DO J=0,NY_S
           IF ((J+INT(RANK/NP_S)*(NY_S+1)).le.NY) THEN
             GY_S(J)=GY(J+INT(RANK/NP_S)*(NY_S+1))
             DY_S(J)=DY(J+INT(RANK/NP_S)*(NY_S+1))
             IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY_S(',J,') = ',GY_S(J)
           END IF
         END DO
         DO K=0,NZ_S
           IF ((K+MOD(RANK,NP_S)*(NZ_S+1)).le.NZ) THEN
             GZ_S(K)=GZ(K+MOD(RANK,NP_S)*(NZ_S+1))
             DZ_S(K)=DZ(K+MOD(RANK,NP_S)*(NZ_S+1))
             IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ_S(',K,') = ',GZ_S(K)
           END IF
         END DO


!         IF ((FLAVOR.NE.'Ensemble') .AND. (RANK.eq.0)) WRITE (6,*)
!     &            'Fourier in X'
         DO I=0,NX_TH
           GX_TH(I)=(I*LX)/NX_TH
           DX_TH(I)=LX/NX_TH
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX_TH(',I,') = ',GX_TH(I)
         END DO
!         IF ((FLAVOR.NE.'Ensemble') .AND. (RANK.eq.0)) WRITE (6,*)
!     &            'Fourier in Z'
         DO K=0,NZ_TH
           GZ_TH(K)=(K*LZ)/NZ_TH
           DZ_TH(K)=LZ/NZ_TH
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ_TH(',K,') = ',GZ_TH(K)
         END DO
!         IF ((FLAVOR.NE.'Ensemble') .AND. (RANK.eq.0)) WRITE (6,*)
!     &            'Fourier in Y'
         DO J=0,NY_TH
           GY_TH(J)=(J*LY)/NY_TH
           DY_TH(J)=LY/NY_TH
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY_TH(',J,') = ',GY_TH(J)
         END DO

         DO J=0,NY_S_TH
           IF ((J+INT(RANK/NP_S)*(NY_S_TH+1)).le.NY_TH) THEN
             GY_S_TH(J)=GY_TH(J+INT(RANK/NP_S)*(NY_S_TH+1))
             DY_S_TH(J)=DY_TH(J+INT(RANK/NP_S)*(NY_S_TH+1))
           END IF
         END DO
         DO K=0,NZ_S_TH
           IF ((K+MOD(RANK,NP_S)*(NZ_S_TH+1)).le.NZ_TH) THEN
             GZ_S_TH(K)=GZ_TH(K+MOD(RANK,NP_S)*(NZ_S_TH+1))
             DZ_S_TH(K)=DZ_TH(K+MOD(RANK,NP_S)*(NZ_S_TH+1))
           END IF
         END DO

         RETURN
         END





C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INPUT_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      REAL    VERSION, CURRENT_VERSION

! Read in input parameters specific for channel flow case
      OPEN (11,file='input_per.dat',form='formatted',status='old')
C Read input file.

      CURRENT_VERSION=1.0
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) VERSION
      IF (VERSION .NE. CURRENT_VERSION)
     &         STOP 'Wrong input data format input_chan'
      if (RANK.eq.0) write(*,*) 'VERSION: ',VERSION
      READ(11,*)
      READ(11,*) TIME_AD_METH
      if (RANK.eq.0) write(*,*) 'TIME_AD_METH: ',TIME_AD_METH
      READ(11,*)
      READ(11,*) LES_MODEL_TYPE
      if (RANK.eq.0) write(*,*) 'LES_MODEL_TYPE: ',LES_MODEL_TYPE
      READ(11,*)
      READ(11,*) IC_TYPE, KICK
      if (RANK.eq.0) write(*,*) 'IC_TYPE,KICK: ',IC_TYPE,KICK
      READ(11,*)
      READ(11,*) I_RO_TAU, PHI, GAMMA, G_TAU, BETA
      if (RANK.eq.0) write(*,*) 'I_RO_TAU,PHI,GAMMA,G_TAU,BETA: '
     &            ,I_RO_TAU,PHI,GAMMA,G_TAU,BETA
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) BACKGROUND_GRAD(N)
        if (RANK.eq.0) write(*,*) 'BACKGROUND_GRAD(N): '
     &                        ,BACKGROUND_GRAD(N)
      END DO

      RETURN
      END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE VIS_FLOW_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS_PER(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      CHARACTER*60 FNAME
      LOGICAL FINAL

      integer i,j,k,n, IZ, IY
      real*8 uc, ubulk

      integer iRANK, RANKY_MOV, RANKZ_MOV

      real*8 ubar,vbar,wbar
      real*8 thbar(1:N_TH)

      real*8 urms_b_sum,vrms_b_sum,wrms_b_sum
      real*8 thth_sum(1:N_TH),thvar_sum(1:N_TH),thme_sum(1:N_TH)
      real*8 l_th(1:N_TH),l_th_sum(1:N_TH)

        IF (USE_MPI) THEN
           CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        END IF

        IF (RANK.eq.0) THEN
          write(*,*) 'TIME , DELTA_T = ',TIME,DELTA_T
        END IF

         IF ((FLAVOR.EQ.'Basic').or.(FLAVOR.eq.'CHEMOTAXIS')) THEN

! Note that this routine uses CFi and CFTH for storage, so it should
! only be used between full R-K timesteps

      IF (RANK.eq.0) THEN
        WRITE(6,*) 'Saving flow statistics.'
      END IF

! Store the velocity in Fourier space in CRi
      do j=0,TNKY
        do k=0,TNKZ_S
          do i=0,NKX_S
            CR1(i,k,j)=CU1(i,k,j)
            CR2(i,k,j)=CU2(i,k,j)
            CR3(i,k,j)=CU3(i,k,j)
          end do
        end do
      end do
! Store the scalar in Fourier space
      do n=1,N_TH
      do j=0,TNKY_TH
        do k=0,TNKZ_S_TH
          do i=0,NKX_S_TH
            CRTH(i,k,j,n)=CTH(i,k,j,n)
          end do
        end do
      end do
      end do

! Compute the vertical gradients in fourier space, store in CFi
      do j=0,TNKY
        do k=0,TNKZ_S
          do i=0,NKX_S
            CF1(i,k,j)=CIKY(j)*CU1(i,k,j)
            CF2(i,k,j)=CIKY(j)*CU2(i,k,j)
            CF3(i,k,j)=CIKY(j)*CU3(i,k,j)
          end do
        end do
      end do
! Save the  scalar in fourier space in CFTH
      do n=1,N_TH
        do j=0,TNKY_TH
          do k=0,TNKZ_S_TH
            do i=0,NKX_S_TH
              CFTH(i,k,j,n)=CTH(i,k,j,n)
            end do
          end do
        end do
      end do

      if (RANK.eq.0) then
        write(*,*) 'Calling FFT'
      end if

        IF (USE_MPI) THEN
           CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        END IF

! Now, convert the velocity and vertical gradients to physical space
      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_PHYSICAL(CU1,U1)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CU2,U2)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CU3,U3)
        do n=1,N_TH
          CALL FFT_XZY_MPI_TO_PHYSICAL_TH(CTH(0,0,0,n),TH(0,0,0,n))
        end do
        CALL FFT_XZY_MPI_TO_PHYSICAL(CF1,F1)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CF2,F2)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CF3,F3)
      ELSE
        CALL FFT_XZY_TO_PHYSICAL(CU1,U1)
        CALL FFT_XZY_TO_PHYSICAL(CU2,U2)
        CALL FFT_XZY_TO_PHYSICAL(CU3,U3)
        do n=1,N_TH
          CALL FFT_XZY_TO_PHYSICAL_TH(CTH(0,0,0,n),TH(0,0,0,n))
        end do
        CALL FFT_XZY_TO_PHYSICAL(CF1,F1)
        CALL FFT_XZY_TO_PHYSICAL(CF2,F2)
        CALL FFT_XZY_TO_PHYSICAL(CF3,F3)
      END IF

       if (RANK.eq.0) then
         write(*,*) 'Done FFT'
       end if

!      if (MOVIE.AND.(RANK.eq.0)) then
! Output a 2d slice through the velocity field for animation in matlab
!        open(79,file='line_v.txt',status='unknown',form='formatted')
!        write(*,*) NX,TIME,DT
!        do i=0,NXM
!          write(79,*) U2(i,NZ_S/2,NY_S/2)
!        end do

!      END IF

      IF (MOVIE) THEN
!      IF (USE_MPI) THEN
!        FNAME='movie_xz_'//MPI_IO_NUM//'.txt'
!      ELSE
!        FNAME='movie_xz.txt'
!      END IF
!      open(80,file=FNAME,status='unknown',form='formatted')
!      do k=0,NZ_S
!      do i=0,NXM
!        write(80,*) GX(i),GZ_S(k),U1(i,k,0)
!      end do
!      end do

! ##### NEW #####
      RANKY_MOV = NY_MOV/(NY_S+1)
      RANKZ_MOV = NZ_MOV/(NZ_S+1)
      call mpi_barrier(MPI_COMM_WORLD,ierror)
      if (RANKY==RANKY_MOV) then
        call WriteHDF5_xzplane(NY_MOV)
      end if
      call mpi_barrier(MPI_COMM_WORLD,ierror)
      if (RANKZ==RANKZ_MOV) then
        call WriteHDF5_xyplane(NZ_MOV)
      end if
      call mpi_barrier(MPI_COMM_WORLD,ierror)
      call WriteHDF5_yzplane(NX_MOV)
! ###############

!      IF (USE_MPI) THEN
!        FNAME='movie_xy_'//MPI_IO_NUM//'.txt'
!      ELSE
!        FNAME='movie_xy.txt'
!      END IF
!      open(81,file=FNAME,status='unknown',form='formatted')
!      do j=0,NY_S
!      do i=0,NXM
!        write(81,*) GX(i),GY_S(j),U1(i,0,j)
!      end do
!      end do



      END IF



! First get the number of samples taken so far
!      NSAMPLES=NSAMPLES+1
! Get the mean velocity



! Get the bulk rms value
      urms_b=0.
      vrms_b=0.
      wrms_b=0.
      do i=0,NXM
!        do j=0,min(NY_S,NYM-((NY_S+1)*(RANK/NP_S)-1))
        do j=0,NY_S
          do k=0,NZ_S
          urms_b=urms_b+(U1(i,k,j)-dble(CR1(0,0,0)))**2.d0
          vrms_b=vrms_b+(U2(i,k,j)-dble(CR2(0,0,0)))**2.d0
          wrms_b=wrms_b+(U3(i,k,j)-dble(CR3(0,0,0)))**2.d0
          end do
        end do
      end do
      urms_b=urms_b/dble(NX*NY*NZ)
      vrms_b=vrms_b/dble(NX*NY*NZ)
      wrms_b=wrms_b/dble(NX*NY*NZ)

      IF (USE_MPI) THEN
        CALL MPI_ALLREDUCE(urms_b,urms_b_sum,1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
        CALL MPI_ALLREDUCE(vrms_b,vrms_b_sum,1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
        CALL MPI_ALLREDUCE(wrms_b,wrms_b_sum,1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
      END IF

      IF (RANK.eq.0) THEN
        write(*,*) '<U_rms>: ',urms_b_sum
        write(*,*) '<V_rms>: ',vrms_b_sum
        write(*,*) '<W_rms>: ',wrms_b_sum
      END IF

      ! Write out the mean statistics at each time
      !       IF (RANK.eq.0) THEN
      !         FNAME='mean.txt'
      !         open(40,file=FNAME,form='formatted',status='unknown')
      !         write(40,*) TIME_STEP,TIME,DELTA_T
      !         write(40,401) urms_b_sum,vrms_b_sum,wrms_b_sum
      !       END IF
      !401   format(3(F20.9,' '))

            IF (RANK.eq.0) THEN
              call WriteStatH5('urms',urms_b_sum)
              call WriteStatH5('vrms',vrms_b_sum)
              call WriteStatH5('wrms',wrms_b_sum)
            END IF

!      if (MOVIE) then
! Output a 2d slice through the velocity field for animation in matlab
!      IF (USE_MPI) THEN
!        FNAME='movie_xz_'//MPI_IO_NUM//'.txt'
!      ELSE
!        FNAME='movie_xz.txt'
!      END IF
!        open(79,file=FNAME,status='unknown'
!     &      ,form='formatted')
!        write(*,*) 'NX_S,NZ_S: ',NX_S,NZ_S
!        do i=0,NXM
!        do k=0,NZ_S
!          write(79,*)
!     &           U1(I,K,10)
!     &           sqrt(U1(I,K,10)**2.d0+U3(I,K,10)**2.d0)
!        end do
!        end do

!      IF (USE_MPI) THEN
!        FNAME='movie_xy_'//MPI_IO_NUM//'.txt'
!      ELSE
!        FNAME='movie_xy.txt'
!      END IF
!        open(80,file=FNAME,status='unknown'
!     &        ,form='formatted')
!        do i=0,NXM
!        do j=0,NY_S
!          write(80,*) U1(I,10,J)
!        end do
!        end do

! This file will contain a single plane and is used in conjunction with
! the matlab script 'realtime_movie' to visualize data during
! simulation
!        open (76,file='temp.txt',status='unknown',form='formatted')
!        do K=0,NZ_S
!          write(76,*) gz(k)
!        end do
!        do I=0,NXM
!        do K=0,NZ_S
!          write(76,*) U1(I,K,NY/4)**2.+U3(I,K,NY/4)**2.+U2(I,K,NY/4)**2.
!        end do
!        end do
!        close (76)
!        CALL SYSTEM('mv temp.txt ../post_process/matlab/latest_slice.txt
!     &')
!      end if

! Do over the number of passive scalars
      do n=1,N_TH
!      do j=0,NY_S
!        thrms(j,n)=0.
!      do k=0,NZ_S
!      do i=0,NXM
!        thrms(j,n)=thrms(j,n)+(abs(TH(i,k,j,n)-mean_th(j,n)))**2.
!      end do
!      end do
!        thrms(j,n)=sqrt(thrms(j,n)/(float(NZ)*float(NX)))
!      end do
! Compute the Reynolds stress and mean velocity gradient
!      do j=0,NY_S
!        thv(j,n)=0.
!      do k=0,NZ_S
!      do i=0,NXM
!       thv(j,n)=thv(j,n)+(TH(i,k,j,n)-mean_th(j,n))
!     +    *(U2(i,k,n)-mean_u2(j))
!      end do
!      end do
!      thv(j,n)=thv(j,n)/(float(NZ)*float(NX))
!      end do

! Get the y-derivative of the mean velocity at GYF points
!      do j=2,NY
!        dthdy(j,n)=(CRTH(0,0,j+1,n)-CRTH(0,0,j-1,n))/(2.*DYF(j))
!      end do
!      do j=1,NY-2
!        dthdy(j,n)=(mean_th(j+1,n)-mean_th(j-1,n))/(GY(j+1)-GY(j-1))
!      end do
!      j=0
!       dthdy(j,n)=(mean_th(j+1,n)-mean_th(NYM,n))/(2.d0*(GY(j+1)-GY(j)))
!      j=NYM
!       dthdy(j,n)=(mean_th(0,n)-mean_th(j-1,n))/(2.d0*(GY(j)-GY(j-1)))

! Get the bacterial/nutrient correlation
      thth(n)=0.d0
      thvar(n)=0.d0
      do j=0,NY_S_TH
      do k=0,NZ_S_TH
      do i=0,NXM_TH
        IF (n.eq.1) THEN
          thth(n)=thth(n)+1.d0*TH(i,k,j,1)
     &                 *DX_TH(I)*DY_TH(J)*DZ_TH(K)
        ELSE
          thth(n)=thth(n)+TH(i,k,j,n)*TH(i,k,j,1)
     &                 *DX_TH(I)*DY_TH(J)*DZ_TH(K)
        END IF
        thvar(n)=thvar(n)+TH(i,k,j,n)*TH(i,k,j,n)
     &                   *DX_TH(I)*DY_TH(J)*DZ_TH(K)
        thme(n)=thme(n)+TH(i,k,j,n)
     &                   *DX_TH(I)*DY_TH(J)*DZ_TH(K)
      end do
      end do
      end do
      thth(n)=thth(n)
      thvar(n)=thvar(n)
      thme(n)=thme(n)/(LX*LY*LZ)

      IF (USE_MPI) THEN
        CALL MPI_ALLREDUCE(THTH(n),THTH_sum(n),1,MPI_DOUBLE_PRECISION
     &                    ,MPI_SUM,MPI_COMM_WORLD,IERROR)
        CALL MPI_ALLREDUCE(THVAR(n),THVAR_sum(n),1,MPI_DOUBLE_PRECISION
     &                    ,MPI_SUM,MPI_COMM_WORLD,IERROR)
        CALL MPI_ALLREDUCE(THME(n),THME_sum(n),1,MPI_DOUBLE_PRECISION
     &                    ,MPI_SUM,MPI_COMM_WORLD,IERROR)
      END IF

      IF (RANK.eq.0) THEN
        write(*,*) 'n,thth(n),thvar(n),thme(n): '
     &         ,n,thth_sum(n),thvar_sum(n)
      END IF


!      if (MOVIE) then
! Output a 2d slice through the scalar field for animation in matlab
!        if (n.eq.1) then
! Chose which scalar is to be outputted
!        open(75,file='movie1_xy.txt',status='unknown',form='formatted')
!        do i=0,NXM
!        do j=0,NY_S
!          write(75,*) TH(I,NZ_S/2,J,n)
!        end do
!        end do
!        open(77,file='movie1_xz.txt',status='unknown',form='formatted')
!        do i=0,NXM
!        do k=0,NZ_S
!          write(77,*) TH(I,K,NY/2,n)
!        end do
!        end do
!        end if
!        if (n.eq.2) then
! Chose which scalar is to be outputted
!        open(85,file='movie2_xy.txt',status='unknown',form='formatted')

!        do i=0,NXM
!        do j=0,NY_S
!          write(85,*) TH(I,NZ_S/2,J,n)
!        end do
!        end do
!        open(87,file='movie2_xz.txt',status='unknown',form='formatted')
!        do i=0,NXM
!        do k=0,NZ_S
!          write(87,*) TH(I,K,NY/2,n)
!        end do
!        end do
!        end if
!        if (n.eq.3) then
! Chose which scalar is to be outputted
!        open(95,file='movie3_xy.txt',status='unknown',form='formatted')
!        do i=0,NXM
!        do j=0,NY_S
!          write(95,*) TH(I,NZ_S/2,J,n)
!        end do
!        end do
!        open(97,file='movie3_xz.txt',status='unknown',form='formatted')
!        do i=0,NXM
!        do k=0,NZ_S
!          write(97,*) TH(I,K,NY/2,n)
!        end do
!        end do
!        end if
!      end if


! End do over number of passive scalars, n
      end do

        IF (USE_MPI) THEN
           CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        END IF

! Convert back to Fourier space
      IF (USE_MPI) THEN
        DO N=1,N_TH
          call FFT_XZY_MPI_TO_FOURIER_TH(TH(0,0,0,n),CTH(0,0,0,n))
        END DO
      ELSE
        DO N=1,N_TH
          call FFT_XZY_TO_FOURIER_TH(TH(0,0,0,n),CTH(0,0,0,n))
        END DO
      END IF


      ! Write out the mean statistics at each time
      !       IF (RANK.eq.0) THEN
      !         FNAME='mean_th.txt'
      !         open(41,file=FNAME,form='formatted',status='unknown')
      !         write(41,*) TIME_STEP,TIME,DELTA_T
      !         do n=1,N_TH
      !         write(41,411) thth_sum(n),thvar_sum(n),thme_sum(n)
      !         end do
      !       END IF
      !411   format(4(F30.15,' '))

            IF (RANK.eq.0) THEN
              call WriteStatH5('thth',thth_sum(1))
              call WriteStatH5('thvar',thvar_sum(1))
              call WriteStatH5('thme',thme_sum(1))
            END IF

!      IF (USE_MPI) THEN
!        FNAME='mean_th'//MPI_IO_NUM//'.txt'
!      ELSE
!        FNAME='mean_th.txt'
!      END IF
!      open(41,file=FNAME,form='formatted',status='unknown')
!      write(41,*) TIME_STEP,TIME,DELTA_T,UBULK
!      do n=1,N_TH
!      do j=0,NYM
!        write(41,402) j,GYF(J),mean_th(j,n)
!     &      ,dthdy(j,n),thrms(j,n),thv(j,n),pe_diss(j,n)
!     &      ,thth(n),thvar(n)
!      end do
!      end do

402   format(I3,' ',8(F30.25,' '))


!      if (VERBOSITY.gt.4) then
!      write(*,*) 'Outputting info for gnuplot...'
!      open (unit=10, file="solution")
!      do i=2,NXM
!        do j=2,NY_S
!          write (10,*) i, j, U1(i,0,j)
!        end do
!        write (10,*) ""
!      end do
!      close (10)
!      call system ('gnuplot <gnuplot.in')
!      end if

! Calculate the TKE
      EK=0.5d0*(urms_b_sum**2.d0+vrms_b_sum**2.d0+wrms_b_sum**2.d0)

        IF (USE_MPI) THEN
           CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        END IF

C Convert velocity back to Fourier space
      IF (USE_MPI) THEN
        call fft_xzy_mpi_to_fourier(U1,CU1)
        call fft_xzy_mpi_to_fourier(U2,CU2)
        call fft_xzy_mpi_to_fourier(U3,CU3)
      ELSE
        call fft_xzy_to_fourier(U1,CU1)
        call fft_xzy_to_fourier(U2,CU2)
        call fft_xzy_to_fourier(U3,CU3)
      END IF

      end if

      call tkebudget_per


      RETURN
      END


      subroutine filter_per(n)
C This subroutine applies a filter to the highest wavenumbers
C It should be applied to the scalars in Fourier space
C The filter used is a sharpened raised cosine filter

      INCLUDE 'header'

      integer I,J,K,N

! Variables for horizontal filtering
      real*8 sigma0

      DO i=0,NKX_S_TH
       DO k=0,TNKZ_S_TH
        DO j=0,TNKY_TH
          sigma0=0.5d0*(1.d0+
     &       cos(sqrt((KX_S_TH(i)*LX*1.d0/float(NX_TH))**2.d0
     &            +(KZ_S_TH(k)*LZ*1.d0/float(NZ_TH))**2.d0
     &            +(KY_TH(j)*LY*1.d0/float(NY_TH))**2.d0)))
! Apply a sharpened raised cosine filter
          CTH(i,k,j,n)=CTH(i,k,j,n)*
     &          sigma0**4.d0*(35.d0-84.d0*sigma0
     &        +70.d0*sigma0**2.d0-20.d0*sigma0**3.d0)
        END DO
       END DO
      END DO

       return
       end


      subroutine tkebudget_per
! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps

      INCLUDE 'header'

      real*8 epsilon_mean,epsilon_sum,eta
      real*8 chi_mean,chi_sum

      integer i,j,k

! Compute the scalar dissipation rate, chi=<dth/dx_i dth/dx_i>
      chi_mean=0.
! Store dth/dx in CSTH1
      do j=0,TNKY_TH
        do k=0,TNKZ_S_TH
          do i=0,NKX_S_TH
            CSTH1(i,k,j)=CIKX_S_TH(i)*CTH(i,k,j,1)
          end do
        end do
      end do
! Convert to physical space
      if (USE_MPI) then
        call FFT_XZY_MPI_TO_PHYSICAL_TH(CSTH1,STH1)
      else
        call FFT_XZY_TO_PHYSICAL_TH(CSTH1,STH1)
      end if
      do j=0,NY_S_TH
        do k=0,NZ_S_TH
          do i=0,NXM_TH
            chi_mean=chi_mean+(STH1(i,k,j)**2.0)
          end do
        end do
      end do
! Store dth/dy in CSTH1
      do j=0,TNKY_TH
        do k=0,TNKZ_S_TH
          do i=0,NKX_S_TH
            CSTH1(i,k,j)=CIKY_TH(j)*CTH(i,k,j,1)
          end do
        end do
      end do
! Convert to physical space
      if (USE_MPI) then
        call FFT_XZY_MPI_TO_PHYSICAL_TH(CSTH1,STH1)
      else
        call FFT_XZY_TO_PHYSICAL_TH(CSTH1,STH1)
      end if
      do j=0,NY_S_TH
        do k=0,NZ_S_TH
          do i=0,NXM_TH
            chi_mean=chi_mean+(STH1(i,k,j)**2.0)
          end do
        end do
      end do
! Store dth/dz in CSTH1
      do j=0,TNKY_TH
        do k=0,TNKZ_S_TH
          do i=0,NKX_S_TH
            CSTH1(i,k,j)=CIKZ_S_TH(k)*CTH(i,k,j,1)
          end do
        end do
      end do
! Convert to physical space
      if (USE_MPI) then
        call FFT_XZY_MPI_TO_PHYSICAL_TH(CSTH1,STH1)
      else
        call FFT_XZY_TO_PHYSICAL_TH(CSTH1,STH1)
      end if
      do j=0,NY_S_TH
        do k=0,NZ_S_TH
          do i=0,NXM_TH
            chi_mean=chi_mean+(STH1(i,k,j)**2.0)
          end do
        end do
      end do
! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
      epsilon_mean=0.
! Store du/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKX_S(i)*CU1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      IF (USE_MPI) THEN
        call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      ELSE
        call FFT_XZY_TO_PHYSICAL(CS1,S1)
      END IF
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKX_S(i)*CU2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      IF (USE_MPI) THEN
        call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      ELSE
        call FFT_XZY_TO_PHYSICAL(CS1,S1)
      END IF
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
      end do
      end do
      end do
! Compute du/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKY(j)*CU1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      IF (USE_MPI) THEN
        call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      ELSE
        call FFT_XZY_TO_PHYSICAL(CS1,S1)
      END IF
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
      end do
      end do
      end do
! Store dw/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKX_S(i)*CU3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      IF (USE_MPI) THEN
        call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      ELSE
        call FFT_XZY_TO_PHYSICAL(CS1,S1)
      END IF
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
      end do
      end do
      end do
! Compute du/dz at GYF gridpoints
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKZ_S(k)*CU1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      IF (USE_MPI) THEN
        call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      ELSE
        call FFT_XZY_TO_PHYSICAL(CS1,S1)
      END IF
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
      end do
      end do
      end do
! Compute dv/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKY(j)*CU2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      IF (USE_MPI) THEN
        call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      ELSE
        call FFT_XZY_TO_PHYSICAL(CS1,S1)
      END IF
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
      end do
      end do
      end do
! Compute dw/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKY(j)*CU3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      IF (USE_MPI) THEN
        call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      ELSE
        call FFT_XZY_TO_PHYSICAL(CS1,S1)
      END IF
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
      end do
      end do
      end do
! Store dv/dz in CS1
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKZ_S(k)*CU2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      IF (USE_MPI) THEN
        call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      ELSE
        call FFT_XZY_TO_PHYSICAL(CS1,S1)
      END IF
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
      end do
      end do
      end do
! Store dw/dz in CS1
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKZ_S(k)*CU3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      IF (USE_MPI) THEN
        call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      ELSE
        call FFT_XZY_TO_PHYSICAL(CS1,S1)
      END IF
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+(S1(i,k,j)**2.0)
      end do
      end do
      end do

      if (RANK.eq.0) write(*,*) 'NX,NY,NZ: ',NX,NY,NZ

      epsilon_mean=NU*epsilon_mean/float(NX*NY*NZ)
      chi_mean=chi_mean/float(NX*NY*NZ)

      IF (USE_MPI) THEN
        CALL MPI_ALLREDUCE(epsilon_mean,epsilon_sum,1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
        CALL MPI_ALLREDUCE(chi_mean,chi_sum,1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
      END IF

      eta=(NU**3.d0/epsilon_sum)**(+0.25d0)

      IF (RANK.eq.0) write(*,*) 'epsilon_mean: ',epsilon_sum

      RE_LAMBDA=EK*SQRT(15.d0/(epsilon_sum*NU))
      IF (RANK.eq.0) write(*,*) 'RE_LAMBDA: ',RE_LAMBDA

      ! Write out the mean statistics at each time
      !      IF (RANK.eq.0) THEN
      !        open(45,file='tke.txt',form='formatted',status='unknown')
      !        write(45,*) TIME_STEP,TIME,DELTA_T
      !        write(45,401) epsilon_sum,eta,re_lambda
      !      END IF
      !401   format(3(F20.9,' '))

            IF (RANK.eq.0) THEN
              call WriteStatH5('epsilon',epsilon_sum)
              call WriteStatH5('eta',eta)
              call WriteStatH5('re_lambda',re_lambda)
              call WriteStatH5('chi',chi_sum)
            END IF

      return
      end
