C periodic.f, the fully-periodic-box solvers for diablo.           VERSION 0.9
C
C These solvers were written by Tom Bewley and John Taylor.
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      INCLUDE 'header'
      
      REAL*8 alpha
      integer i,j,k

C Define variables for the Geophysical case
      PI=4.D0*ATAN(1.D0)
      PHI=2.D0*PI*PHI/360.D0
      GAMMA=2.D0*PI*GAMMA/360.D0
      C_SIN=COS(PHI)*SIN(GAMMA)/SIN(PHI)
      C_COS=COS(PHI)*COS(GAMMA)/SIN(PHI)

      if (F_TYPE.eq.3) then
        do j=0,TNKY
          do k=0,TNKZ_S
            do i=0,NKX_S
              call RANDOM_NUMBER(alpha)
              f_phase(i,k,j)=alpha*2.d0*pi
            end do
          end do
        end do
      end if

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
      TEMP1=NU*H_BAR(RK_STEP)/2.d0
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

      ! Set the CFx variables to zero here so they can be used as working
      ! variables in USER_RHS_PER_FOURIER
      CF1(:,:,:)=0.d0
      CF2(:,:,:)=0.d0
      CF3(:,:,:)=0.d0

      ! Add the forcing terms CFx, CFTH using the subroutine in user_rhs.f
      ! Note that CUx, CTH are stored in Fourier space here
      call USER_RHS_PER_FOURIER


C For Geophysical applications, add the Coriolis term
C Note, on an f-plane under the traditional approximation, C_SIN=C_COS=0
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
           CF1(I,K,J)=CF1(I,K,J)+I_RO_TAU*
     &                  (C_SIN*CU2(I,K,J)+CU3(I,K,J))
           CF2(I,K,J)=CF2(I,K,J)+I_RO_TAU*
     &                  (-1.d0*C_COS*CU3(I,K,J)+C_SIN*CU1(I,K,J))
           CF3(I,K,J)=CF3(I,K,J)+I_RO_TAU*
     &                  (C_COS*CU2(I,K,J)-CU1(I,K,J))
!           CF1(I,K,J)=I_RO_TAU*CU3(I,K,J)
!           CF2(I,K,J)=0.d0
!           CF3(I,K,J)=-1.d0*I_RO_TAU*CU1(I,K,J)
          END DO
        END DO
      END DO


C Transform the scalar concentration to physical space
      DO N=1,N_TH
        CALL FFT_XZY_MPI_TO_PHYSICAL_TH(CTH(0,0,0,N),TH(0,0,0,N))
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
        CALL FFT_XZY_MPI_TO_FOURIER_TH(STH1,CSTH1)
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
        CALL FFT_XZY_MPI_TO_FOURIER_TH(STH1,CSTH1)
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
        CALL FFT_XZY_MPI_TO_FOURIER_TH(STH1,CSTH1)
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
      CALL FFT_XZY_MPI_TO_PHYSICAL(CU1,U1)
      CALL FFT_XZY_MPI_TO_PHYSICAL(CU2,U2)
      CALL FFT_XZY_MPI_TO_PHYSICAL(CU3,U3)

C Compute the nonlinear terms for the momentum equations
        DO J=0,NY_S
          DO K=0,NZ_S
            DO I=0,NXM
              S1(I,K,J)=U1(I,K,J)*U1(I,K,J)
            END DO
          END DO
        END DO
        CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
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
        CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
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
        CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
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
        CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
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
        CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
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
        CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
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
        CALL FFT_XZY_MPI_TO_FOURIER_TH(TH(0,0,0,N),CTH(0,0,0,N))
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
        CALL COURANT_MPI
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

      CALL FFT_XZY_MPI_TO_PHYSICAL(CR1,R1)
      CALL FFT_XZY_MPI_TO_PHYSICAL(CR2,R2)

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
      CALL FFT_XZY_MPI_TO_PHYSICAL(CF1,F1)
      CALL FFT_XZY_MPI_TO_PHYSICAL(CF2,F2)
      CALL FFT_XZY_MPI_TO_PHYSICAL(CF3,F3)
      CALL FFT_XZY_MPI_TO_PHYSICAL(CR1,R1)
      CALL FFT_XZY_MPI_TO_PHYSICAL(CR2,R2)
      CALL FFT_XZY_MPI_TO_PHYSICAL(CR3,R3)
      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            P(I,K,J)=2*(P(I,K,J)+F1(I,K,J)*F2(I,K,J)
     *                          +F3(I,K,J)*R1(I,K,J)
     *                          +R2(I,K,J)*R3(I,K,J))
          END DO
        END DO
      END DO
      CALL FFT_XZY_MPI_TO_FOURIER(P,CP)

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
      REAL*8 K0,K_MAG, a0, theta0, phi0
      CHARACTER*60 FNAME
      REAL*8 b3,f03,E3,Sigma3,A3,alpha,h3,beta3,N03


C For an initial vortex, define the location of the centerline
      REAL*8 XC(0:NY+1),ZC(0:NY+1)

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

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
          CU1=0.d0
          CU2=0.d0
          CU3=0.d0
          if (RANK.eq.0) then
            CU1(1,0,1)=0.1/sqrt(2.0)/2.0
!            CU1(0,0,1)=0.5
            CU2(1,0,1)=-0.1/sqrt(2.0)/2.0
            CU1(0,2,0) = 0.5
          end if
        ELSE IF (IC_TYPE.eq.3) THEN
        ! Initialize with a GM spectrum of internal waves
          b3=1300./1e-2/sqrt(RI_TAU(1)/NU)
          f03=7.3e-5*1e2*sqrt(RI_TAU(1))
          N03=5.2e-3*1e2*sqrt(RI_TAU(1))
          E3=6.3e-5
          Sigma3=0.468
          A3=sqrt(RI_TAU(1)*b3*E3*f03*sqrt(RI_TAU(1)-f03**2))
     &          /sqrt(2*PI*Sigma3)
          h3=1/sqrt(2.)
          do j=0,TNKY
            do k=0,TNKZ_S
              do i=0,NKX_S
                if ( (KX2_S(i)+KZ2_S(k)+KY2(j).le.100.) .and.
     &              (KY(j).ne.0)) then
                  if (KX2_S(i)+KZ2_S(k).eq.0) then ! Shear Flow component
                    CS1(i,k,j)=sqrt(2*b3*E3*RI_TAU(1)/PI/Sigma3/h3**2)
     &    *(KY2(j)+9*RI_TAU(1)*pi**2/b3**2/N03**2)**(-0.5)*(atan(h3*
     &        sqrt(RI_TAU(1)-f03**2)/f03/sqrt(h3**2+KY2(j))))**(0.5)
                    call RANDOM_NUMBER(beta3)
                    call RANDOM_NUMBER(alpha)
                    CU1(i,k,j)=sqrt(beta3)*cexp(cmplx(0,2.*pi*alpha))
     &                      *CS1(i,k,j)
                    call RANDOM_NUMBER(alpha)
                    CU3(i,k,j)=sqrt(1-beta3)*cexp(cmplx(0,2.*pi*alpha))
     &                      *CS1(i,k,j)
                  else
                    call RANDOM_NUMBER(alpha)
                    alpha=2.*pi*alpha ! Random phase of each wave
                    CS1(i,k,j)=A3*cexp(cmplx(0,alpha))*KY(j)/
     &         (sqrt(KX2_S(i)+KZ2_S(k))*
     &          sqrt(KX2_S(i)+KZ2_S(k)+KY2(j))
     &         *(RI_TAU(1)*(KX2_S(i)+KZ2_S(k))+f03**2*KY2(j))
     &         *(KY2(j)+pi**2*9*RI_TAU(1)/N03**2/b3**2))**(0.5)
                    CU1(i,k,j)=CS1(i,k,j)*KY(j)*KX_S(i)/sqrt(
     &            (KX2_S(i)+KZ2_S(k)+KY2(j))*(KX2_S(i)+KZ2_S(k)))
                    CU2(i,k,j)=CS1(i,k,j)*sqrt(KX2_S(i)+KZ2_S(k))/sqrt(
     &            KX2_S(i)+KZ2_S(k)+KY2(j))
                    CU3(i,k,j)=CS1(i,k,j)*KY(j)*KZ_S(k)/sqrt(
     &            (KX2_S(i)+KZ2_S(k)+KY2(j))*(KX2_S(i)+KZ2_S(k)))
                    CTH(i,k,j,1)=CS1(i,k,j)*CI/sqrt(RI_TAU(1))
                  end if
                end if
              end do
            end do
          end do
        else if (IC_TYPE.eq.4) then
          ! Initialise with large shear + monochromatic IGW
          a0 = 0.1  ! steepness
          i = 2     ! KX
          j = 4     ! KY
          CU1=0.d0
          CU2=0.d0
          CU3=0.d0
          if (RANK.eq.0) then
            CU1(0,0,1)=-CI ! Vertical shear at wavenumber 1
            CU1(i,0,j)=CI*a0/2*sqrt(RI_TAU(1))/
     &                    sqrt(real(i**2+j**2))
            CU2(i,0,j)=-CI*real(i)/real(j)*a0/2*sqrt(RI_TAU(1))/
     &                    sqrt(real(i**2+j**2))
          end if
        ELSE
          WRITE(*,*) 'Warning, Undefined ICs in periodic.f'
        END IF

        if (IC_TYPE.lt.2) then
          CALL FFT_XZY_MPI_TO_FOURIER(U1,CU1)
          CALL FFT_XZY_MPI_TO_FOURIER(U3,CU3)
          CALL FFT_XZY_MPI_TO_FOURIER(U2,CU2)
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

      ELSE
        write(6,*) 'Unknown flavor, flow-field not created'
      end if

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

      if (RANK.eq.0) write(*,*) 'calling rem_div...'
      CALL REM_DIV_PER

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K,N
      REAL*8 a0,theta0,phi0

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
          CTH(:,:,:,N)=0.d0
            if ((RANK.eq.0) .and. (N.eq.1)) then
              CTH(1,0,1,N)=CI*0.1/2.0
            end if
          END IF
        END DO
      else if (IC_TYPE.eq.4) then
        a0 = 0.1  ! steepness
        i = 2     ! KX
        j = 4     ! KY
        if (CREATE_NEW_TH(1) .and. RANK.eq.0) then
          CTH(i,0,j,1)=-a0/2/real(j)
        end if
      ELSE
        WRITE(*,*) 'UNKNOWN IC_TYPE IN CREATE_TH_PER'
      END IF

! Transfer to Fourier space
      if (IC_TYPE.ne.2 .and. IC_TYPE.ne.4) then
        DO N=1,N_TH
          CALL FFT_XZY_MPI_TO_FOURIER_TH(TH(0,0,0,N),CTH(0,0,0,N))
        END DO
      end if

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_GRID_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K

      DO I=0,NX
        GX(I)=(I*LX)/NX
        DX(I)=LX/NX
        IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
      END DO

      DO K=0,NZ
        GZ(K)=(K*LZ)/NZ
        DZ(K)=LZ/NZ
        IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ(',K,') = ',GZ(K)
      END DO

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


         DO I=0,NX_TH
           GX_TH(I)=(I*LX)/NX_TH
           DX_TH(I)=LX/NX_TH
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX_TH(',I,') = ',GX_TH(I)
         END DO

         DO K=0,NZ_TH
           GZ_TH(K)=(K*LZ)/NZ_TH
           DZ_TH(K)=LZ/NZ_TH
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ_TH(',K,') = ',GZ_TH(K)
         END DO

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
      INTEGER N

! Read in input parameters specific for periodic flow case
      OPEN (11,file='input_per.dat',form='formatted',status='old')
C Read input file.

      CURRENT_VERSION=2.0
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) VERSION
      IF (VERSION .NE. CURRENT_VERSION)
     &         STOP 'Wrong input data format input_chan'
      if (RANK.eq.0) write(*,*) 'VERSION: ',VERSION
      READ(11,*)
      READ(11,*) IC_TYPE, KICK
      if (RANK.eq.0) write(*,*) 'IC_TYPE,KICK: ',IC_TYPE,KICK
      READ(11,*)
      READ(11,*) I_RO_TAU, PHI, GAMMA, G_TAU, BETA
      if (RANK.eq.0) write(*,*) 'I_RO_TAU,PHI,GAMMA,G_TAU,BETA: '
     &            ,I_RO_TAU,PHI,GAMMA,G_TAU,BETA
      READ(11,*)
      READ(11,*) F_TYPE, FORCE_SHEAR, target_Reb
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) BACKGROUND_GRAD(N)
        if (RANK.eq.0) write(*,*) 'BACKGROUND_GRAD(N): '
     &                        ,BACKGROUND_GRAD(N)
      END DO

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

      real*8 U1rms,U2rms,U3rms
      real*8 U1rms_sum,U2rms_sum,U3rms_sum
      real*8 THrms_sum(1:N_TH)
      real*8 l_th(1:N_TH),l_th_sum(1:N_TH),thflux_sum(1:N_TH)

      real*8 U1me(0:NY_S), U3me(0:NY_S), TH1me(0:NY_S)
      real*8 U1U2(0:NY_S), U3U2(0:NY_S), THU2(0:NY_S_TH,1:N_TH)
      real*8 U1U1(0:NY_S), U2U2(0:NY_S), U3U3(0:NY_S)
      real*8 U1rms_h(0:NY_S), U2rms_h(0:NY_S), U3rms_h(0:NY_S)
      real*8 THTH(0:NY_S_TH,1:N_TH), THTH_h(0:NY_S_TH,1:N_TH)
      real*8 U1U2_sum(0:NY_S), U3U2_sum(0:NY_S)
      real*8 THU2_sum(0:NY_S_TH,1:N_TH)
      real*8 E_L(0:TNKY), E_S(0:TNKY), spectrum(0:TNKY), E(0:TNKY)
      character(10) :: gname

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

      IF (RANK.eq.0) THEN
        write(*,*) 'TIME , DELTA_T = ',TIME,DELTA_T
      END IF

      IF ((FLAVOR.EQ.'Basic').or.(FLAVOR.eq.'CHEMOTAXIS')) THEN

! Note that this routine uses CRi and CRTH for storage, so it should
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

      if (RANK.eq.0) then
        write(*,*) 'Calling FFT'
      end if

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

! Now, convert the velocity and vertical gradients to physical space
      CALL FFT_XZY_MPI_TO_PHYSICAL(CU1,U1)
      CALL FFT_XZY_MPI_TO_PHYSICAL(CU2,U2)
      CALL FFT_XZY_MPI_TO_PHYSICAL(CU3,U3)
      do n=1,N_TH
        CALL FFT_XZY_MPI_TO_PHYSICAL_TH(CTH(0,0,0,n),TH(0,0,0,n))
      end do

      if (RANK.eq.0) then
        write(*,*) 'Done FFT'
      end if

      IF (MOVIE) THEN
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
      END IF

! Get the horizontal means of U1, U3 & TH
      call horiz_mean(CR1,U1me)
      call horiz_mean(CR3,U3me)
      call horiz_mean(CRTH,TH1me)
      call mpi_barrier(MPI_COMM_WORLD,ierror)
      if (rank.eq.0) write(*,*) 'Horizontal means computed'

      if (RANKZ.eq.0) then
      gname='U1me'
      call WriteMeanH5(gname,U1me)
      gname='U3me'
      call WriteMeanH5(gname,U3me)
      gname='THme'
      call WriteMeanH5(gname,TH1me)
      end if
      call mpi_barrier(MPI_COMM_WORLD,ierror)

      if (rank.eq.0) write(*,*) 'Horizontal means written to file'
     
! Get the bulk rms values and write perturbation velocities to CRx
      U1rms=0.d0
      U2rms=0.d0
      U3rms=0.d0
      U1U2=0.d0
      U3U2=0.d0
      U1U1=0.d0
      U2U2=0.d0
      U3U3=0.d0

      do i=0,NXM
        do j=0,NY_S
          do k=0,NZ_S
            S1(i,k,j)=U1(i,k,j)-U1me(j)
            U1rms=U1rms+S1(i,k,j)**2.d0
            U1U1(j)=U1U1(j)+S1(i,k,j)**2.d0
            U1U2(j)=U1U2(j)+S1(i,k,j)*U2(i,k,j)
          end do
        end do
      end do
      call fft_xzy_mpi_to_fourier(S1,CR1)
      do i=0,NXM
        do j=0,NY_S
          do k=0,NZ_S
            S1(i,k,j)=U2(i,k,j)-dble(CR2(0,0,0))
            U2rms=U2rms+S1(i,k,j)**2.d0
            U2U2(j)=U2U2(j)+S1(i,k,j)**2.d0
            S1(i,k,j)=U3(i,k,j)-U3me(j)
            U3rms=U3rms+S1(i,k,j)**2.d0
            U3U3(j)=U3U3(j)+S1(i,k,j)**2.d0
            U3U2(j)=U3U2(j)+S1(i,k,j)*U2(i,k,j)
          end do
        end do
      end do
      call fft_xzy_mpi_to_fourier(S1,CR3)

      U1rms=U1rms/dble(NX*NY*NZ)
      U2rms=U2rms/dble(NX*NY*NZ)
      U3rms=U3rms/dble(NX*NY*NZ)
      U1U2=U1U2/dble(NX*NZ)
      U3U2=U3U2/dble(NX*NZ)
      U1U1=U1U1/dble(NX*NZ)
      U2U2=U2U2/dble(NX*NZ)
      U3U3=U3U3/dble(NX*NZ)

      CALL MPI_ALLREDUCE(U1rms,U1rms_sum,1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(U2rms,U2rms_sum,1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(U3rms,U3rms_sum,1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(U1U2,U1U2_sum,NY_S+1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,IERROR)
      CALL MPI_ALLREDUCE(U3U2,U3U2_sum,NY_S+1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,IERROR)
      CALL MPI_ALLREDUCE(U1U1,U1rms_h,NY_S+1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,IERROR)
      CALL MPI_ALLREDUCE(U2U2,U2rms_h,NY_S+1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,IERROR)
      CALL MPI_ALLREDUCE(U3U3,U3rms_h,NY_S+1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,IERROR)

      IF (RANK.eq.0) THEN
        U1rms_sum=sqrt(U1rms_sum)
        U2rms_sum=sqrt(U2rms_sum)
        U3rms_sum=sqrt(U3rms_sum)

        write(*,*) '<U1rms>: ',U1rms_sum
        write(*,*) '<U2rms>: ',U2rms_sum
        write(*,*) '<U3rms>: ',U3rms_sum
        
        gname='U1rms'
        call WriteStatH5(gname,U1rms_sum)
        gname='U2rms'
        call WriteStatH5(gname,U2rms_sum)
        gname='U3rms'
        call WriteStatH5(gname,U3rms_sum)
      END IF

      if (RANKZ.eq.0) then
        gname='U1U2'
        call WriteMeanH5(gname,U1U2_sum)
        gname='U3U2'
        call WriteMeanH5(gname,U3U2_sum)
        U1rms_h=sqrt(U1rms_h)
        U2rms_h=sqrt(U2rms_h)
        U3rms_h=sqrt(U3rms_h)
        gname='U1rms'
        call WriteMeanH5(gname,U1rms_h)
        gname='U2rms'
        call WriteMeanH5(gname,U2rms_h)
        gname='U3rms'
        call WriteMeanH5(gname,U3rms_h)
      end if
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

! Do over the number of passive scalars
      do n=1,N_TH

! Get the bacterial/nutrient correlation
      THrms(n)=0.d0
      thflux(n)=0.d0
      THU2(:,n)=0.d0
      do j=0,NY_S_TH
      do k=0,NZ_S_TH
      do i=0,NXM_TH
        THrms(n)=THrms(n)+TH(i,k,j,n)*TH(i,k,j,n)
        thflux(n)=thflux(n)+TH(i,k,j,n)*U2(i,k,j)
        THU2(j,n)=THU2(j,n)+TH(i,k,j,n)*U2(i,k,j)
        THTH(j,n)=THTH(j,n)+(TH(i,k,j,n)-TH1me(j))**2
      end do
      end do
      end do
      THrms(n)=THrms(n)/dble(NX_TH*NY_TH*NZ_TH)
      thflux(n)=thflux(n)/dble(NX_TH*NY_TH*NZ_TH)
      THU2(:,n)=THU2(:,n)/dble(NX_TH*NZ_TH)
      THTH(:,n)=THTH(:,n)/dble(NX_TH*NZ_TH)

      CALL MPI_ALLREDUCE(THrms(n),THrms_sum(n),1,MPI_DOUBLE_PRECISION
     &                    ,MPI_SUM,MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(THFLUX(n),THFLUX_sum(n),1
     &        ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(THU2,THU2_sum,NY_S+1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,IERROR)
      CALL MPI_ALLREDUCE(THTH,THTH_h,NY_S+1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,IERROR)

      IF (RANK.eq.0) THEN
        THrms_sum=sqrt(THrms_sum)

        write(*,*) 'n,THrms(n),THflux(n): '
     &         ,n,THrms_sum(n),THflux_sum(n)

        gname='THrms'
        call WriteStatH5(gname,THrms_sum(n))
        gname='THflux'
        call WriteStatH5(gname,THflux_sum(n))
      END IF

      if (RANKZ.eq.0) then
      gname='THflux'
      call WriteMeanH5(gname,THU2_sum(:,n))
      THTH_h=sqrt(THTH_h)
      gname='THrms'
      call WriteMeanH5(gname,THTH_h(:,n))
      end if
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

! End do over number of passive scalars, n
      end do

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

! Convert back to Fourier space
      DO N=1,N_TH
        call FFT_XZY_MPI_TO_FOURIER_TH(TH(0,0,0,n),CTH(0,0,0,n))
      END DO

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

C Convert velocity back to Fourier space
        call fft_xzy_mpi_to_fourier(U1,CU1)
        call fft_xzy_mpi_to_fourier(U2,CU2)
        call fft_xzy_mpi_to_fourier(U3,CU3)

      end if


      call spectra_per(CU1,E)
      spectrum=0.d0
      call MPI_ALLREDUCE(E,spectrum,TNKY+1,MPI_DOUBLE_PRECISION
     &                    ,MPI_SUM,MPI_COMM_WORLD,IERROR)
      gname='U1'
      if (RANK.EQ.0) call WriteSpectrumH5(gname,spectrum)

      call spectra_per(CU2,E)
      spectrum=0.d0
      call MPI_ALLREDUCE(E,spectrum,TNKY+1,MPI_DOUBLE_PRECISION
     &                    ,MPI_SUM,MPI_COMM_WORLD,IERROR)
      gname='U2'
      if (RANK.EQ.0) call WriteSpectrumH5(gname,spectrum)
      
      call spectra_per(CU3,E)
      spectrum=0.d0
      call MPI_ALLREDUCE(E,spectrum,TNKY+1,MPI_DOUBLE_PRECISION
     &                    ,MPI_SUM,MPI_COMM_WORLD,IERROR)
      gname='U3'
      if (RANK.EQ.0) call WriteSpectrumH5(gname,spectrum)
      
      call spectra_per(CTH(:,:,:,1),E)
      spectrum=0.d0
      call MPI_ALLREDUCE(E,spectrum,TNKY+1,MPI_DOUBLE_PRECISION
     &                    ,MPI_SUM,MPI_COMM_WORLD,IERROR)
      gname='TH1'
      if (RANK.EQ.0) call WriteSpectrumH5(gname,spectrum)
      
      call tkebudget_per


      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      subroutine tkebudget_per
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps

      INCLUDE 'header'

      real*8 epsilon_mean,epsilon_sum,eta
      real*8 chi_mean,chi_sum
      real*8 epsilon_h(0:NY_S), chi_h(0:NY_S)
      real*8 epsilon_h_sum(0:NY_S), chi_h_sum(0:NY_S)
      character(10) :: gname

      integer i,j,k

! Compute the scalar dissipation rate, chi=<dth/dx_i dth/dx_i>
      chi_mean=0.d0
      chi_h(:)=0.d0
! Store dth/dx in CSTH1
      do j=0,TNKY_TH
        do k=0,TNKZ_S_TH
          do i=0,NKX_S_TH
            CSTH1(i,k,j)=CIKX_S_TH(i)*CTH(i,k,j,1)
          end do
        end do
      end do
! Convert to physical space
      call FFT_XZY_MPI_TO_PHYSICAL_TH(CSTH1,STH1)
      do j=0,NY_S_TH
        do k=0,NZ_S_TH
          do i=0,NXM_TH
            chi_mean=chi_mean+(STH1(i,k,j)**2.0)
            chi_h(j)=chi_h(j)+(STH1(i,k,j)**2.0)
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
      call FFT_XZY_MPI_TO_PHYSICAL_TH(CSTH1,STH1)
      do j=0,NY_S_TH
        do k=0,NZ_S_TH
          do i=0,NXM_TH
            chi_mean=chi_mean+(STH1(i,k,j)**2.0)
            chi_h(j)=chi_h(j)+(STH1(i,k,j)**2.0)
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
      call FFT_XZY_MPI_TO_PHYSICAL_TH(CSTH1,STH1)
      do j=0,NY_S_TH
        do k=0,NZ_S_TH
          do i=0,NXM_TH
            chi_mean=chi_mean+(STH1(i,k,j)**2.0)
            chi_h(j)=chi_h(j)+(STH1(i,k,j)**2.0)
          end do
        end do
      end do
! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
      epsilon_mean=0.d0
      epsilon_h=0.d0
! Store du/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKX_S(i)*CR1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+(S1(i,k,j)**2.0)
        epsilon_h(j)=epsilon_h(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKX_S(i)*CR2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
        epsilon_h(j)=epsilon_h(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKY(j)*CR1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
        epsilon_h(j)=epsilon_h(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dw/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKX_S(i)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
        epsilon_h(j)=epsilon_h(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dz at GYF gridpoints
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKZ_S(k)*CR1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
        epsilon_h(j)=epsilon_h(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute dv/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKY(j)*CR2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
        epsilon_h(j)=epsilon_h(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute dw/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKY(j)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
        epsilon_h(j)=epsilon_h(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dz in CS1
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKZ_S(k)*CR2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+S1(i,k,j)**2.0
        epsilon_h(j)=epsilon_h(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dw/dz in CS1
      do j=0,TNKY
      do k=0,TNKZ_S
      do i=0,NKX_S
        CS1(i,k,j)=CIKZ_S(k)*CR3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        epsilon_mean=epsilon_mean+(S1(i,k,j)**2.0)
        epsilon_h(j)=epsilon_h(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do

      if (RANK.eq.0) write(*,*) 'NX,NY,NZ: ',NX,NY,NZ

      epsilon_mean=NU*epsilon_mean/float(NX*NY*NZ)
      chi_mean=NU/Pr(1)*chi_mean/float(NX*NY*NZ)

      CALL MPI_ALLREDUCE(epsilon_mean,epsilon_sum,1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(chi_mean,chi_sum,1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)

      CALL MPI_ALLREDUCE(epsilon_h,epsilon_h_sum,NY_S+1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,IERROR)
      CALL MPI_ALLREDUCE(chi_h,chi_h_sum,NY_S+1
     &              ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,IERROR)

      epsilon_h_sum=NU*epsilon_h_sum/float(NX*NZ)
      chi_h_sum=NU/Pr(1)*chi_h_sum/float(NX*NZ)

      IF (RANK.eq.0) write(*,*) 'epsilon_mean: ',epsilon_sum

            IF (RANK.eq.0) THEN
              gname='epsilon'
              call WriteStatH5(gname,epsilon_sum)
              gname='chi'
              call WriteStatH5(gname,chi_sum)
            END IF

      if (RANKZ.eq.0) then
      gname='epsilon'
      call WriteMeanH5(gname,epsilon_h_sum)
      gname='chi'
      call WriteMeanH5(gname,chi_h_sum)
      end if
      call mpi_barrier(MPI_COMM_WORLD,ierror)

      return
      end

C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
      subroutine spectra_per(CU,E)
C----*|--.---------.---------.---------.---------.---------.---------.-|-----|
! WORK IN PROGRESS SPECTRUM CALCULATOR
      INCLUDE 'header'
      
      complex*16 CU(0:NX_S/2,0:NZ_S,0:NY+1)
      real*8 E(0:TNKY)
      integer i,j,k

      E=0.d0
      CS1=0.5*CU*conjg(CU)
      do j=0,TNKY
        do k=0,TNKZ_S
          do i=0,NKX_S
            if ((KX_S(I).EQ.0) .AND. ((KZ_S(K).LT.0) .OR. 
     &            ((KZ_S(K).EQ.0) .AND. (KY(J).LT.0)))) then
            else
              if ((i+RANKZ*(NKX_S+1).LE.NKX) .AND.
     &            (k+RANKY*(TNKZ_S+1).LE.TNKZ)) then
                if ((RANK.EQ.0).AND.
     &            (i.eq.0) .AND. (j.eq.0) .AND. (k.eq.0)) then
                  E(j)=E(j)+real(CS1(i,k,j))
                else
                  E(j)=E(j)+2.d0*real(CS1(i,k,j))
                end if
              end if
            end if
            
          end do
        end do
      end do

      end subroutine spectra_per
