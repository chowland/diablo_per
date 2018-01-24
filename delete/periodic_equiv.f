
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
!      EPSILON_TARGET=((1.d0/(DX(1)*0.25d0))**4.d0)
!     &              *(NU**3.d0)*(PR(1))**(-2.d0)
      EPSILON_TARGET=((1.d0/DX(1))**4.d0)*(NU**3.d0)*(100)**(-2.d0)
      write(*,*) 'EPSILON_TARGET: ',EPSILON_TARGET
!      END IF		   

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_PER_1
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Main time-stepping algorithm for the fully periodic case
C INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE "mpif.h"

      INCLUDE 'header'
      INCLUDE 'header_mpi'

      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6
      INTEGER I,J,K,N

C Start with data local in the X-direction

C Compute the RHS in Fourier Space, CRi.

C First, define the constants used for time-stepping
      TEMP1=NU*H_BAR(RK_STEP)/2.0
      TEMP2=BETA_BAR(RK_STEP)*H_BAR(RK_STEP)
      TEMP3=ZETA_BAR(RK_STEP)*H_BAR(RK_STEP)
      TEMP4=H_BAR(RK_STEP)
      TEMP5=0.0



      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
C Start with the explicit part of the Crank-Nicolson viscous term and
C  the pressure gradient treated with Explicit Euler:
            TEMP5=1-TEMP1*(KX2_S(I)+KY2(J)+KZ2_S(K))

            CR1(I,K,J)=TEMP5*CU1(I,K,J)-TEMP4*(CIKX_S(I)*CP(I,K,J))
            CR2(I,K,J)=TEMP5*CU2(I,K,J)-TEMP4*(CIKY(J)*CP(I,K,J))
            CR3(I,K,J)=TEMP5*CU3(I,K,J)-TEMP4*(CIKZ_S(K)*CP(I,K,J))
C For each scalar, start with the explict part of the Crank-Nicolson
C diffusive term for each scalar
            DO N=1,N_TH
              TEMP6=1-(TEMP1/PR(N))*(KX2_S(I)+KY2(J)+KZ2_S(K))
     &               -REACTION(N)*TEMP1
              CRTH(I,K,J,N)=TEMP6*CTH(I,K,J,N)
            END DO
          END DO
        END DO
        IF (RK_STEP .GT. 1) THEN
          DO K=0,TNKZ_S
            DO I=0,NKX_S
C Add the term: ZETA_BAR(RK_STEP)*R(U(RK_STEP-1))
              CR1(I,K,J)=CR1(I,K,J)+TEMP3*CF1(I,K,J)
              CR2(I,K,J)=CR2(I,K,J)+TEMP3*CF2(I,K,J)
              CR3(I,K,J)=CR3(I,K,J)+TEMP3*CF3(I,K,J)
C Do the same for each scalar:
              DO N=1,N_TH
                CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP3*CFTH(I,K,J,N)
              END DO
            END DO
          END DO
        END IF
      END DO
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
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CFTH(I,K,J,N)=-CU2(I,K,J)
          END DO
        END DO
      END DO 
      ELSE
! Otherwise don't
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CFTH(I,K,J,N)=0.D0
          END DO
        END DO
      END DO 
      END IF
      END DO



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


C Here, for the Chemotaxis simulations,
C   transform only the bacterial concentration to physical space
C   leave the nutrient concentration (TH(:,:,:,1)) in Fourier space
      DO N=2,N_TH 
        IF (USE_MPI) THEN
          CALL FFT_XZY_MPI_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N))
        ELSE
          CALL FFT_XZY_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N))
        END IF
      END DO

C Calculate the Chemotaxis terms
C Note that it is assumed that the nutrient concentration is in scalar #1
      IF (FLAVOR.eq.'CHEMOTAXIS') THEN
      IF (N_TH.ge.1) THEN
C If we are using the SPEED-LIMITING form:
C First calculate tanh(k|grad(c)|)/|grad(c)|
C Note that we can use F1 and F2 as temporary storage variables here
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CF1(I,K,J)=CIKX_S(I)*CTH(I,K,J,1)
          END DO
        END DO
      END DO
      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_PHYSICAL(CF1,F1)
      ELSE
        CALL FFT_XZY_TO_PHYSICAL(CF1,F1)
      END IF
      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            F2(I,K,J)=F1(I,K,J)**2.d0
          END DO
        END DO
      END DO
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CF1(I,K,J)=CIKY(J)*CTH(I,K,J,1)
          END DO
        END DO
      END DO
      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_PHYSICAL(CF1,F1)
      ELSE
        CALL FFT_XZY_TO_PHYSICAL(CF1,F1)
      END IF
      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            F2(I,K,J)=F2(I,K,j)+F1(I,K,J)**2.d0
          END DO
        END DO
      END DO
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CF1(I,K,J)=CIKZ_S(K)*CTH(I,K,J,1)
          END DO
        END DO
      END DO
      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_PHYSICAL(CF1,F1)
      ELSE
        CALL FFT_XZY_TO_PHYSICAL(CF1,F1)
      END IF
! Note, here the constant k depends on the particular value of CHI
      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            F2(I,K,J)=F2(I,K,j)+F1(I,K,J)**2.d0
            F2(I,K,J)=sqrt(F2(I,K,J))
            F1(I,K,J)=V_S(2)*TANH(V_S(2)*F2(I,K,J)/CHI(2))/F2(I,K,J)
          END DO
        END DO
      END DO
C F1 now contains the saturation term
C Do for each bacterial species
      DO N=2,N_TH
      IF (CHI(N).NE.0.d0) THEN
C Here, calculate each of the three chemotaxis terms
C First, CHI*d/dx(B*dC/dx):
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CS1(I,K,J)=CIKX_S(I)*CTH(I,K,J,1)
          END DO
        END DO
      END DO
      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      ELSE
        CALL FFT_XZY_TO_PHYSICAL(CS1,S1)
      END IF
      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            S1(I,K,J)=S1(I,K,J)*TH(I,K,J,N)*F1(I,K,J)
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
            CFTH(I,K,J,N)=CFTH(I,K,J,N)-CIKX_S(I)*CS1(I,K,J)
          END DO
        END DO
      END DO
 
C Now, calculate the CHI*d/dy(B*dC/dy) term:
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CS1(I,K,J)=CIKY(J)*CTH(I,K,J,1)
          END DO
        END DO
      END DO
      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      ELSE
        CALL FFT_XZY_TO_PHYSICAL(CS1,S1)
      END IF
      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            S1(I,K,J)=S1(I,K,J)*TH(I,K,J,N)*F1(I,K,J)
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
            CFTH(I,K,J,N)=CFTH(I,K,J,N)-CIKY(J)*CS1(I,K,J)
          END DO
        END DO
      END DO
C Finally, calculate the CHI*d/dz(B*dC/dz) term:
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CS1(I,K,J)=CIKZ_S(K)*CTH(I,K,J,1)
          END DO 
        END DO
      END DO
      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_PHYSICAL(CS1,S1)
      ELSE
        CALL FFT_XZY_TO_PHYSICAL(CS1,S1)
      END IF
      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            S1(I,K,J)=S1(I,K,J)*TH(I,K,J,N)*F1(I,K,J)
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
            CFTH(I,K,J,N)=CFTH(I,K,J,N)-CIKZ_S(K)*CS1(I,K,J)
          END DO
        END DO
      END DO
C End if CHI(N).NE.0
      END IF
C End loop over bacterial species
      END DO
C End if N_TH>=1
      END IF
C End if CHEMOTAXIS
      END IF

C Now that we are done with the Chemotaxis terms, we need to transform the
C nutrient concentration to physical space for calculating the nonlinear
C advection terms

      IF (N_TH.gt.0) THEN
        IF (USE_MPI) THEN
          CALL FFT_XZY_MPI_TO_PHYSICAL(CTH(0,0,0,1),TH(0,0,0,1))
        ELSE
          CALL FFT_XZY_TO_PHYSICAL(CTH(0,0,0,1),TH(0,0,0,1))
        END IF
      END IF

C Now the velocity and scalar fields are in physical space and local in y

C Add the nutrient uptake terms
      IF (FLAVOR.eq.'CHEMOTAXIS') THEN
        DO N=2,N_TH
          DO j=0,NY_S
            DO k=0,NZ_S
              DO i=0,NXM
                S1(I,K,J)=C0(N)*TH(I,K,J,1)*TH(I,K,J,N)
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
                CFTH(I,K,J,1)=CFTH(I,K,J,1)-CS1(I,K,J)
              END DO
            END DO
          END DO
        END DO
      END IF   


C Compute the nonlinear terms for the passive scalar equation
C Do this before the nonlinear momentum terms to use Fi as a working
C  array before using it for the momentum equation.


      DO N=1,N_TH
        DO J=0,NY_S
          DO K=0,NZ_S
            DO I=0,NXM
              F1(I,K,J)=U1(I,K,J)*TH(I,K,J,N)
              F2(I,K,J)=U2(I,K,J)*TH(I,K,J,N)
              F3(I,K,J)=U3(I,K,J)*TH(I,K,J,N)
            END DO
          END DO
        END DO
        IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_FOURIER(F1,CF1)
        CALL FFT_XZY_MPI_TO_FOURIER(F2,CF2)
        CALL FFT_XZY_MPI_TO_FOURIER(F3,CF3)
        ELSE
        CALL FFT_XZY_TO_FOURIER(F1,CF1)
        CALL FFT_XZY_TO_FOURIER(F2,CF2)
        CALL FFT_XZY_TO_FOURIER(F3,CF3)
        END IF
        DO J=0,TNKY
          DO K=0,TNKZ_S
            DO I=0,NKX_S
              CFTH(I,K,J,N)=CFTH(I,K,J,N)-CIKX_S(I)*CF1(I,K,J)
     &                      -CIKY(J)*CF2(I,K,J)
     &                      -CIKZ_S(K)*CF3(I,K,J)
     
C Add R-K terms for the TH equation to the RHS
              CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP2*CFTH(I,K,J,N)
            END DO
          END DO
        END DO
      END DO
C The RHS vector for the TH equation is now ready



C Compute the nonlinear terms for the momentum equations
      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            F1(I,K,J)=U1(I,K,J)*U1(I,K,J)
            F2(I,K,J)=U1(I,K,J)*U2(I,K,J)
            F3(I,K,J)=U1(I,K,J)*U3(I,K,J)
            S1(I,K,J)=U2(I,K,J)*U2(I,K,J)
          END DO
        END DO
      END DO


      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_FOURIER(F1,CF1)
        CALL FFT_XZY_MPI_TO_FOURIER(F2,CF2)
        CALL FFT_XZY_MPI_TO_FOURIER(F3,CF3)
        CALL FFT_XZY_MPI_TO_FOURIER(S1,CS1)
      ELSE
        CALL FFT_XZY_TO_FOURIER(F1,CF1)
        CALL FFT_XZY_TO_FOURIER(F2,CF2)
        CALL FFT_XZY_TO_FOURIER(F3,CF3)
        CALL FFT_XZY_TO_FOURIER(S1,CS1)
      END IF


C Here we start constructing the R-K terms in CFi
C Note, that the order of the following operations are important
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            CF1(I,K,J)=-CIKX_S(I)*CF1(I,K,J)
     *                 -CIKY(J)*CF2(I,K,J)
     *                 -CIKZ_S(K)*CF3(I,K,J)
            CF2(I,K,J)=-CIKX_S(I)*CF2(I,K,J)
     *                 -CIKY(J)*CS1(I,K,J)
            CF3(I,K,J)=-CIKX_S(I)*CF3(I,K,J)
          END DO
        END DO
       END DO
C At this point, F1,F2,F3 contain some of the nonlinear terms


C Compute the remaining nonlinear terms
C We cannot use F1,F2,F3 for storage, so use S1 as the working variable
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
C Fisrt, convert back to Fourier space
       IF (USE_MPI) THEN
         CALL FFT_XZY_MPI_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N))
       ELSE
         CALL FFT_XZY_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N))
       END IF
       DO J=0,TNKY
         DO K=0,TNKZ_S
           DO I=0,NKX_S
             CF2(I,K,J)=CF2(I,K,J)+RI_TAU(N)*CTH(I,K,J,N)
           END DO
         END DO
       END DO
       END DO

C NEEDS TO BE PARALLELIZED!!!!!
C Add some forcing to the system to keep the Batchelor scale fixed
      EK_local=0.d0
      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
              EK_local=EK_local
     &            +U1(I,K,J)**2.d0+U2(I,K,J)**2.d0+U3(I,K,J)**2.d0
          END DO
        END DO
      END DO
      IF (USE_MPI) THEN
        CALL MPI_ALLREDUCE(EK_local,EK,1,MPI_DOUBLE_PRECISION,MPI_SUM
     &                     ,MPI_COMM_WORLD,IERROR)
      END IF 
C Note, that each cell has the same volume, so we can just average over all points
      EK=EK/dble(NX*NY*NZ)

      EPSILON_TARGET=15.d0*EK**2.d0/(35.d0**2.d0*NU)

! Scale EK by an amount to compensate for dissipation from 2/3 de-aliasing:
      EK=0.8d0*EK

      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            S1(I,K,J)=(EPSILON_TARGET/EK)*U1(I,K,J)
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
            CF1(I,K,J)=CF1(I,K,J)+CS1(I,K,J)
          END DO
        END DO
      END DO
      DO J=0,NY_S
       DO K=0,NZ_S
          DO I=0,NXM
            S1(I,K,J)=(EPSILON_TARGET/EK)*U2(I,K,J)
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
            CF2(I,K,J)=CF2(I,K,J)+CS1(I,K,J)
          END DO
        END DO
      END DO
      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            S1(I,K,J)=(EPSILON_TARGET/EK)*U3(I,K,J)
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
            CF3(I,K,J)=CF3(I,K,J)+CS1(I,K,J)
          END DO
        END DO
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

C Now solve the implicit system for the intermediate field.
C (In the fully-periodic case, this is easy!)
      DO J=0,TNKY
        DO K=0,TNKZ_S
          DO I=0,NKX_S
            TEMP5=1.d0+TEMP1*(KX2_S(I)+KY2(J)+KZ2_S(K))
            CU1(I,K,J)=CR1(I,K,J)/TEMP5
            CU2(I,K,J)=CR2(I,K,J)/TEMP5
            CU3(I,K,J)=CR3(I,K,J)/TEMP5
            DO N=1,N_TH
              TEMP6=1+(TEMP1/PR(N))*(KX2_S(I)+KY2(J)+KZ2_S(K)) 
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
      INCLUDE "mpif.h"
      INCLUDE 'header'
      INCLUDE 'header_mpi'
      INTEGER I,J,K
      REAL*8  TEMP5

C Compute phi, store in the variable CR1.
C Note the coefficient H_BAR is absorbed into phi.
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
      CP(0,0,0)=0.0

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'mpif.h'
      INCLUDE 'header'
      INCLUDE 'header_mpi'
      INTEGER I, J, K
      REAL*8 RNUM1,RNUM2,RNUM3
      REAL*8 K0
      CHARACTER*35 FNAME

C For an initial vortex, define the location of the centerline
      REAL*8 XC(0:NY+1),ZC(0:NY+1)

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

      IF ((FLAVOR .EQ. 'Basic').OR.(FLAVOR.eq.'CHEMOTAXIS')) THEN

      WRITE(6,*) 'Creating new flow from scratch.'

C Initialize random number generator
      CALL RANDOM_SEED

      IF (IC_TYPE.eq.0) THEN
C Initizlize the flow using a Taylor-Green vortex
C Nondimensionalize with U0 and 1/kappa

      DO J=0,NY_S
        DO K=0,NZ_S
          DO I=0,NXM
            IF (((J+INT(RANK/NP_S)*(NY_S+1)).le.NY)
     &          .and.((K+MOD(RANK,NP_S)*(NZ_S+1)).le.NZ)) then
! Add a random phase
            U1(I,K,J)=cos(2*pi*(GY_S(J))/LY)
     &               *cos(2*pi*(GX(I))/LX)
     &               *SIN(2*pi*(GZ_S(K))/LZ)
!            U1(I,K,J)=cos(8.*pi*GY_S(J)/LY)*sin(2*pi*GX(I)/LX)+1.d0
            U2(I,K,J)=0.d0
!            U3(I,K,J)=0.d0
            U3(I,K,J)=-cos(2*pi*(GY_S(J))/LY)
     &               *sin(2*pi*(GX(I))/LX)
     &               *COS(2*pi*(GZ_S(K))/LZ) 
            else
             U1(i,k,j)=0.d0
             U2(i,k,j)=0.d0
             U3(i,k,j)=0.d0
            end if
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
      ELSE
        WRITE(*,*) 'Warning, Undefined Initial conditions in periodic.f'
      END IF


!      IF (USE_MPI) THEN
!        FNAME='movie_xy_'//MPI_IO_NUM//'.txt'
!      ELSE
!        FNAME='movie_xy.txt'
!      END IF
!        open(80,file=FNAME,status='unknown'
!     &        ,form='formatted')
!        do i=0,NX_S
!        do j=0,NY
!          write(80,*) U1(I,10,J)
!        end do
!        end do

      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_FOURIER(U1,CU1)
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

        CALL FFT_XZY_MPI_TO_FOURIER(U3,CU3)
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

        CALL FFT_XZY_MPI_TO_FOURIER(U2,CU2)
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
       
      ELSE
        CALL FFT_XZY_TO_FOURIER(U1,CU1)
        CALL FFT_XZY_TO_FOURIER(U2,CU2)
        CALL FFT_XZY_TO_FOURIER(U3,CU3)
      END IF

!      DO J=0,NY_S
!        DO K=0,NZ_S
!          DO I=0,NKX
!            CALL RANDOM_NUMBER(RNUM1)
!            CALL RANDOM_NUMBER(RNUM2)
!            CALL RANDOM_NUMBER(RNUM3)
!            CU1(I,K,J)=CU1(I,K,J)+(RNUM1-0.5)*KICK
!            CU2(I,K,J)=CU2(I,K,J)+(RNUM2-0.5)*KICK
!            CU3(I,K,J)=CU3(I,K,J)+(RNUM3-0.5)*KICK
!          end do
!        end do
!      end do
	
	ELSE IF (FLAVOR .EQ. 'Ensemble') THEN
	
	CALL CREATE_FLOW_ENSEM
	
	ELSE
	write(6,*) 'Unknown flavor, flow-field not created'
	
	endif

! get the initial energy in low wavenumbers
!      IF (USE_MPI) THEN
!        write(*,*) 'RANK,CU1: ',RANK,CU1(5,5,5)
!        CALL FFT_XZY_MPI_TO_PHYSICAL(CU1,U1)
!        write(*,*) 'RANK,U1: ',RANK,U1(5,5,5)
!        CALL FFT_XZY_MPI_TO_PHYSICAL(CU2,U2)
!        CALL FFT_XZY_MPI_TO_PHYSICAL(CU3,U3)
!      ELSE
!        CALL FFT_XZY_TO_PHYSICAL(CU1,U1)
!        CALL FFT_XZY_TO_PHYSICAL(CU2,U2)
!        CALL FFT_XZY_TO_PHYSICAL(CU3,U3)
!      END IF
!      EK0=0.d0
!      DO J=0,NYM
!        DO K=0,NZ_S
!          DO I=0,NX_S
!              EK0=EK0+U1(I,K,J)**2.d0+U2(I,K,J)**2.d0+U3(I,K,J)**2.d0
!          END DO
!        END DO
!      END DO
!      write(*,*) 'EK0: ',EK0
!      IF (N_TH.gt.0) THEN
!      EPSILON_TARGET=((1.d0/DX(1))**4.d0)*(NU**3.d0)*(PR(1))**(-2.d0)
!      EPSILON_TARGET=((1.d0/DX(1))**4.d0)*(NU**3.d0)*(100.d0)**(-2.d0)
!      write(*,*) 'EPSILON_TARGET: ',EPSILON_TARGET
!      write(*,*) 'TARGET KOLMOGOROV SCALE: ',
!     &         (NU**3.d0/epsilon_target)**(0.25d0)
!      END IF
!      IF (USE_MPI) THEN
!        CALL FFT_XZY_MPI_TO_FOURIER(U1,CU1)
!        CALL FFT_XZY_MPI_TO_FOURIER(U2,CU2)
!        CALL FFT_XZY_MPI_TO_FOURIER(U3,CU3)
!      ELSE
!        CALL FFT_XZY_TO_FOURIER(U1,CU1)
!        CALL FFT_XZY_TO_FOURIER(U2,CU2)
!        CALL FFT_XZY_TO_FOURIER(U3,CU3)
!      END IF
 
      CALL SAVE_STATS_PER(.FALSE.)


      write(*,*) 'calling rem_div'
      CALL REM_DIV_PER
       write(*,*) 'In CREATE 3 RANK; ',RANK,CU1(10,10,10)


!      CALL POISSON_P_PER

      CALL SAVE_STATS_PER(.FALSE.)
      write(*,*) 'DONE CREATE_FLOW '

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K,N
      real*8 lp

C Note, Since stratification is not permitted in the periodic flow field
C Any background stratification must be added to the governing equations

C Set the patch size in terms of the domain size (lp/L)
!      lp=0.0108695d0
      lp=0.1d0

      DO N=1,N_TH
      IF (CREATE_NEW_TH(N)) THEN
        DO J=0,NY_S
          DO K=0,NZ_S
            DO I=0,NXM
              IF (N.eq.1) THEN
C Nutrient concentration will be centered in the domain
         TH(I,K,J,N)=EXP(-((GX(I)-LX/2)**2.d0/(2.d0*(lp/2.d0)**2.d0))
     &                   -((GY_S(J)-LY/2)**2.d0/(2.d0*(lp/2.d0)**2.d0))
     &                   -((GZ_S(K)-LZ/2)**2.d0/(2.d0*(lp/2.d0)**2.d0)))
C Concentration centered at specific location
!         TH(I,K,J,N)=EXP(-((GX(I)-0.473)**2.d0/(2.d0*(lp/2.d0)**2.d0))
!     &                   -((GY_S(J)-0.316)**2.d0/(2.d0*(lp/2.d0)**2.d0))
!     &                 -((GZ_S(K)-0.402)**2.d0/(2.d0*(lp/2.d0)**2.d0)) )
              ELSE
! Uniform initial distribution
        TH(I,K,J,N)=1.d0

              END IF
            END DO
          END DO
        END DO

       IF (USE_MPI) THEN
         CALL FFT_XZY_MPI_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N)) 
       ELSE
         CALL FFT_XZY_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N)) 
       END IF
      
       END IF
       END DO

       RETURN
       END
      

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_GRID_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'mpif.h'
      INCLUDE 'header' 
      INCLUDE 'header_mpi'
      INTEGER I,J,K

         IF (FLAVOR.NE.'Ensemble') WRITE (6,*) 'Fourier in X'
         DO I=0,NX
           GX(I)=(I*LX)/NX
           DX(I)=LX/NX
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
         END DO
         IF (FLAVOR.NE.'Ensemble') WRITE (6,*) 'Fourier in Z'
         DO K=0,NZ
           GZ(K)=(K*LZ)/NZ
           DZ(K)=LZ/NZ
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ(',K,') = ',GZ(K)
         END DO
         IF (FLAVOR.NE.'Ensemble') WRITE (6,*) 'Fourier in Y'
         DO J=0,NY
           GY(J)=(J*LY)/NY
           DY(J)=LY/NY
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY(',J,') = ',GY(J)
         END DO

         DO J=0,NY_S
           IF ((J+INT(RANK/NP_S)*(NY_S+1)).le.NY) THEN
             GY_S(J)=GY(J+INT(RANK/NP_S)*(NY_S+1))
             DY_S(J)=DY(J+INT(RANK/NP_S)*(NY_S+1))
           END IF
         END DO
         DO K=0,NZ_S
           IF ((K+MOD(RANK,NP_S)*(NZ_S+1)).le.NZ) THEN
             GZ_S(K)=GZ(K+MOD(RANK,NP_S)*(NZ_S+1))
             DZ_S(K)=DZ(K+MOD(RANK,NP_S)*(NZ_S+1))
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
      READ(11,*)
      READ(11,*) TIME_AD_METH
      READ(11,*)
      READ(11,*) LES_MODEL_TYPE
      READ(11,*)
      READ(11,*) IC_TYPE, KICK
      READ(11,*)
      READ(11,*) TRUCK, GUSTS
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) BACKGROUND_GRAD(N)
        READ(11,*)
        READ(11,*) CHI(N), C0(N), V_S(N)
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
      INCLUDE 'mpif.h'
      INCLUDE 'header'
      INCLUDE 'header_mpi'
      CHARACTER*35 FNAME
      LOGICAL FINAL
 
      integer jmax,kmax

      integer i,j,k,n
      real*8 uc, ubulk

      real*8 urms_b_sum,vrms_b_sum,wrms_b_sum
      real*8 thth_sum(1:N_TH),thvar_sum(1:N_TH)

        IF (USE_MPI) THEN
           CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        END IF


         IF ((FLAVOR.EQ.'Basic').or.(FLAVOR.eq.'CHEMOTAXIS')) THEN

! Note that this routine uses CFi and CFTH for storage, so it should
! only be used between full R-K timesteps

      IF (RANK.eq.0) THEN
        WRITE(6,*) 'Saving flow statistics.'
      END IF

      if (FINAL) then
! We are done with the simulation
! Write out statistics to a file
        open(20,file='stats.txt',form='formatted',status='unknown')
        do j=0,NYM
          write(20,201) j,GYF(j),UBAR(j),VBAR(j),WBAR(j)
        end do
201     format(I3,',',F16.9,',',F16.9,',',F16.9,',',F16.9)
        do n=1,N_TH
        do j=0,NYM
          write(20,202) j,GYF(j),THBAR(j,n)
        end do
        end do
202     format(I3,',',F16.9,',',F16.9)

        else
! We are in the middle of a run, compile statistics

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
        do j=0,TNKY
          do k=0,TNKZ_S
            do i=0,NKX_S
              CFTH(i,k,j,n)=CTH(i,k,j,n)
            end do
          end do
        end do
      end do

! Now, convert the velocity and vertical gradients to physical space
      IF (USE_MPI) THEN
        CALL FFT_XZY_MPI_TO_PHYSICAL(CU1,U1)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CU2,U2)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CU3,U3)
        do n=1,N_TH
          CALL FFT_XZY_MPI_TO_PHYSICAL(CTH(0,0,0,n),TH(0,0,0,n))
        end do
        CALL FFT_XZY_MPI_TO_PHYSICAL(CF1,F1)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CF2,F2)
        CALL FFT_XZY_MPI_TO_PHYSICAL(CF3,F3)
      ELSE
        CALL FFT_XZY_TO_PHYSICAL(CU1,U1)
        CALL FFT_XZY_TO_PHYSICAL(CU2,U2)
        CALL FFT_XZY_TO_PHYSICAL(CU3,U3)
        do n=1,N_TH
          CALL FFT_XZY_TO_PHYSICAL(CTH(0,0,0,n),TH(0,0,0,n))
        end do
        CALL FFT_XZY_TO_PHYSICAL(CF1,F1)
        CALL FFT_XZY_TO_PHYSICAL(CF2,F2)
        CALL FFT_XZY_TO_PHYSICAL(CF3,F3)
      END IF

! First get the number of samples taken so far
      NSAMPLES=NSAMPLES+1
! Get the mean velocity

      IF (RANK.eq.0) then
        jmax=NY_S
        kmax=NZ_S
      else if (RANK.eq.1) then
        jmax=NY_S
        kmax=30
      else if (RANK.eq.2) then
        jmax=30
        kmax=NZ_S
      else
        jmax=30
        kmax=30  
      end if
       


! Get the bulk rms value
      urms_b=0.
      vrms_b=0.
      wrms_b=0.
      do i=0,NXM
!        do j=0,min(NY_S,NYM-((NY_S+1)*(RANK/NP_S)-1))
        do j=0,NY_S
          do k=0,NZ_S
          urms_b=urms_b+U1(i,k,j)**2.d0
          vrms_b=vrms_b+U2(i,k,j)**2.d0
          wrms_b=wrms_b+U3(i,k,j)**2.d0
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
!      IF (USE_MPI) THEN
!        FNAME='mean'//MPI_IO_NUM//'.txt'
!      ELSE
!        FNAME='mean.txt'
!      END IF
!      open(40,file=FNAME,form='formatted',status='unknown')
!      write(40,*) TIME_STEP,TIME,DELTA_T,UBULK
!      do j=0,NY_S
!        write(40,401) j,GY(J),mean_u1(j)
!     +      ,mean_u2(j)
!     +      ,mean_u3(j),urms(j),vrms(j),wrms(j)
!     +      ,uv(j),uw(j),wv(j),dudy(j),dwdy(j),mean_p(j),shear(j)
!      end do

401   format(I3,' ',14(F20.9,' '))

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
      do j=0,NY_S
      do k=0,NZ_S
      do i=0,NXM
        thth(n)=thth(n)+TH(i,k,j,n)*TH(i,k,j,1)*DX(I)*DY(J)*DZ(K)
        thvar(n)=thvar(n)+TH(i,k,j,n)*TH(i,k,j,n)*DX(I)*DY(J)*DZ(K)
      end do
      end do
      end do
      thth(n)=thth(n)
      thvar(n)=thvar(n)

      IF (USE_MPI) THEN
        CALL MPI_ALLREDUCE(THTH(n),THTH_sum(n),1,MPI_DOUBLE_PRECISION
     &                    ,MPI_SUM,MPI_COMM_WORLD,IERROR)
        CALL MPI_ALLREDUCE(THVAR(n),THVAR_sum(n),1,MPI_DOUBLE_PRECISION
     &                    ,MPI_SUM,MPI_COMM_WORLD,IERROR)
      END IF

      IF (RANK.eq.0) THEN
        write(*,*) 'n,thth(n),thvar(n): ',n,thth_sum(n),thvar_sum(n)
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


! Convert back to Fourier space
      IF (USE_MPI) THEN
        do n=1,N_TH
          call fft_xzy_mpi_to_fourier(TH(0,0,0,n),CTH(0,0,0,n))
        end do
      ELSE
        do n=1,N_TH
          call fft_xzy_to_fourier(TH(0,0,0,n),CTH(0,0,0,n))
        end do
      END IF


! Write out the mean statistics at each time
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

!      call tkebudget_per
	
	END IF

      RETURN
      END


      subroutine filter_per
C This subroutine applies a filter to the highest wavenumbers
C It should be applied to the scalars in Fourier space
C The filter used is a sharpened raised cosine filter

      include 'header'

      integer I,J,K,N

! Variables for horizontal filtering
      real*8 sigma0

C Set the filtering constants for the all directions
      DO N=1,N_TH
      DO i=0,NKX_S
       DO k=0,TNKZ_S
        DO j=0,TNKY
          sigma0=0.5d0*(1.d0+
     &       cos(sqrt((KX(i)*LX*1.d0/float(NX))**2.d0
     &            +(KZ(k)*LZ*1.d0/float(NZ))**2.d0
     &            +(KY(j)*LY*1.d0/float(NY))**2.d0)))
! Apply a sharpened raised cosine filter
          CTH(i,k,j,n)=CTH(i,k,j,n)*
     &          sigma0**4.d0*(35.d0-84.d0*sigma0
     &        +70.d0*sigma0**2.d0-20.d0*sigma0**3.d0)
        END DO
       END DO
      END DO
      END DO

       return
       end


      subroutine tkebudget_per
! NOte, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps

      include 'header'

      integer i,j,k



! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
      do j=0,NYM
        epsilon(j)=0.
      end do
! Store du/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CU1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CU2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKY(j)*CU1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CF1,F1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dvdx*dudy
        epsilon(j)=epsilon(j)+(S1(i,k,j)*F1(i,k,j))
      end do
      end do
      end do
! Store dw/dx in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKX(i)*CU3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute du/dz at GYF gridpoints, note remove mean
! Store du/dz in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKZ(k)*CU1(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CF1,F1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dudz*dwdx
        epsilon(j)=epsilon(j)+S1(i,k,j)*F1(i,k,j)
      end do
      end do
      end do
! Compute dv/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKY(j)*CU2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Compute dw/dy at GYF gridpoints, note remove mean
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKY(j)*CU3(i,k,j)
      end do
      end do
      end do
! Convert to physical space 
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(S1(i,k,j)**2.0)
      end do
      end do
      end do
! Store dv/dz in CF1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CF1(i,k,j)=CIKZ(k)*CU2(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CF1,F1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+0.5*(F1(i,k,j)**2.0)
! Cross term dvdz*dwdy
        epsilon(j)=epsilon(j)+S1(i,k,j)*F1(i,k,j)
      end do
      end do
      end do
! Store dw/dz in CS1
      do j=0,TNKY
      do k=0,TNKZ
      do i=0,NKX
        CS1(i,k,j)=CIKZ(k)*CU3(i,k,j)
      end do
      end do
      end do
! Convert to physical space
      call fft_xzy_to_physical(CS1,S1)
      do j=0,NYM
      do k=0,NZM
      do i=0,NXM
        epsilon(j)=epsilon(j)+(S1(i,k,j)**2.0)
      end do
      end do
      end do
      do j=0,NYM
        epsilon(j)=epsilon(j)/float(NX*NZ)
      end do


! Write out the bulk rms velocity
      write(*,*) '<U_rms>2: ',urms_b


! Write out the mean statistics at each time
      open(45,file='tke.txt',form='formatted',status='unknown')
      write(45,*) TIME_STEP,TIME,DELTA_T
      do j=0,NYM
        write(45,401) j,GY_S(J),epsilon(j)
      end do
401   format(I3,' ',2(F20.9,' '))


! Get Kolmogorov wavelength
      epsilon_mean=NU*SUM(epsilon(0:NYM))/dble(NY) 

      k_eta=2.d0*PI*(NU**3.d0/epsilon_mean)**(-0.25d0)

      write(*,*) 'Kolmogorov scale: ',2.d0*PI/k_eta

      return
      end






