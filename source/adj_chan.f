C******************************************************************************|
C AD_CHAN.f, the channel-flow adjoint solvers for diablo.           VERSION 0.9
C This solver was written by Chris Colburn (fall 2007).
C******************************************************************************|
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INIT_ADJ_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

C The RK scheme needs to be initialized at CUiS == 0 before this
C algorithm is used.  Make sure that this is properly implemented in the
C main code of diablo.

C Initialize any constants here
      PI=4.D0*ATAN(1.D0)

      RETURN
      END
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE RK_ADJ_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Main time-stepping algorithm for the channel-adjoint case.
C This algorithm uses Crank-Nicolson for all terms involving vertical
C derivatives (viscous and nonlinear) and 3rd order Runge-Kutta for the
C rest of the terms
C INPUTS  (in Fourier space): CUiS, CPS, CUi 
C         and (if k>1) CFi at (k-1)  (for i=1,2,3)
C OUTPUTS (in Fourier space): CUiS, CPS, and (if k<3) CFi at (k)
C Each RK step, there are 14 FFT calls. 11 storage variables are used.     
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header_batch'

      INTEGER I,J,K,N      
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, UBULK

C Define the constants that are used in the time-stepping
C For reference, see Numerical Renaissance
      TEMP1=NU * H_BAR(RK_STEP) / 2.0
      TEMP2=H_BAR(RK_STEP) / 2.0
      TEMP3=ZETA_BAR(RK_STEP) * H_BAR(RK_STEP)
      TEMP4=H_BAR(RK_STEP)
      TEMP5=BETA_BAR(RK_STEP) * H_BAR(RK_STEP)


C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!  INCLUDE PREVIOUS ADJOINT STATES  !!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Store the old velocity in the RHS vector
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CU1(I,K,J)
            CR3(I,K,J)=CU3(I,K,J)
          END DO
        END DO
      END DO
      DO J=2,NY 
        DO K=0,TNKZ
          DO I=0,NKX
            CR2(I,K,J)=CU2(I,K,J)
          END DO
        END DO
      END DO


C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!  PREVIOUS NON-LINEAR TERMS  !!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Add the R-K term from the rk-1 step 
      IF (RK_STEP .GT. 1) THEN
        DO J=2,NYM
          DO K=0,TNKZ
            DO I=0,NKX
              CR1(I,K,J)=CR1(I,K,J)+TEMP3*CF1(I,K,J)
              CR3(I,K,J)=CR3(I,K,J)+TEMP3*CF3(I,K,J)
            END DO
          END DO
        END DO
        DO J=2,NY
          DO K=0,TNKZ
            DO I=0,NKX
              CR2(I,K,J)=CR2(I,K,J)+TEMP3*CF2(I,K,J)
            END DO
          END DO
        END DO
      END IF


C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!  PRESSURE TERMS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Take the y-derivative of the pressure at GY points in Fourier space
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CS1(I,K,J)=(CP(I,K,J) - CP(I,K,J-1)) / DY(J)
          END DO
        END DO
      END DO

C Add the pressure gradient to the RHS as explicit Euler
      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CR1(I,K,J)-TEMP4*(CIKX(I)*CP(I,K,J))
            CR3(I,K,J)=CR3(I,K,J)-TEMP4*(CIKZ(K)*CP(I,K,J))
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CR2(I,K,J)=CR2(I,K,J)-TEMP4*CS1(I,K,J)
          END DO
        END DO
      END DO


C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!  LINEAR TERMS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C COMPUTE d2__d(x1)2 && d2__d(x3)2

      DO J=2,NYM
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CR1(I,K,J) + TEMP1*((KX2(I)+KZ2(K))*CU1(I,K,J))
            CR3(I,K,J)=CR3(I,K,J) + TEMP1*((KX2(I)+KZ2(K))*CU3(I,K,J))
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CR2(I,K,J)=CR2(I,K,J) + TEMP1*(KX2(I)+KZ2(K))*CU2(I,K,J)
          END DO
        END DO
      END DO

C Compute the vertical viscous term in physical space and add to RHS
C This is the explicit part of the Crank-Nicolson term
	DO J=2,NYM
	  DO K=0,TNKZ
	    DO I=0,NKX
          CR1(I,K,J)=CR1(I,K,J)+TEMP1*
     &      (  ((CU1(I,K,J+1) - CU1(I,K,J)) / DY(J+1)  
     &         -(CU1(I,K,J)   - CU1(I,K,J-1)) / DY(J)) /DYF(J)  )
          CR3(I,K,J)=CR3(I,K,J)+TEMP1*
     &      (  ((CU3(I,K,J+1) - CU3(I,K,J)) / DY(J+1) 
     &         -(CU3(I,K,J)   - CU3(I,K,J-1)) / DY(J)) /DYF(J)  )
          END DO
        END DO
      END DO

	DO J=2,NYM
	  DO K=0,TNKZ
	      DO I=0,NKX
            CR2(I,K,J)=CR2(I,K,J)+TEMP1*
     &        (  ((CU2(I,K,J+1) - CU2(I,K,J))  / DYF(J) 
     &           -(CU2(I,K,J)   - CU2(I,K,J-1))/ DYF(J-1))/DY(J)  )
          END DO
        END DO
      END DO


C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!  NON-LINEAR TERMS  !!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Convert flow velocities to physical space, because non-linear products
C in the adjoint equation need to be calculated in physical space.
      call fft_xz_to_physical(CU1F1,U1F1,0,NY+1)
      call fft_xz_to_physical(CU2F1,U2F1,0,NY+1)
      call fft_xz_to_physical(CU3F1,U3F1,0,NY+1)
	
C U1: -U3F*(dU1dX3 + dU3dX1)   &&   U3: -U1F*(dU1dX3 + dU3dX1)

	DO J=2,NYM
	  DO K=0,TNKZ
	    DO I=0,NKX
	    CS1(I,K,J)=CIKZ(K)*CU1(I,K,J)+CIKX(I)*CU3(I,K,J)
	    END DO
	  END DO
	END DO
	
      call fft_xz_to_physical(CS1,S1,0,NY+1)
	
	DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
	      F1(I,K,J) = -U3F1(I,K,J)*S1(I,K,J)
	      F3(I,K,J) = -U1F1(I,K,J)*S1(I,K,J)
	    END DO
	  END DO
	END DO	

C U1: -2*U1F*dU1dX1

	DO J=2,NYM
	  DO K=0,TNKZ
	    DO I=0,NKX
	      CS1(I,K,J) = CIKX(I)*CU1(I,K,J)
	    END DO
	  END DO
	END DO
	
      call fft_xz_to_physical(CS1,S1,0,NY+1)
	
	DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
	      F1(I,K,J) = F1(I,K,J) - 2.*U1F1(I,K,J)*S1(I,K,J)
	    END DO
	  END DO
      END DO
	
	
C U1: -U2F*(dU1dX2 + dU2dX1) ??? Interpolations of U1 & U2

	DO J=2,NYM
	  DO K=0,TNKZ
	    DO I=0,NKX
	      CS1(I,K,J) = ( (((CU1(I,K,J+1)+CU1(I,K,J))*0.5 
     &	  	         - (CU1(I,K,J)+CU1(I,K,J-1))*0.5)/DYF(J))
     &			   + (CIKX(I)*(CU2(I,K,J+1)+CU2(I,K,J))/2.0)  )
	    END DO
	  END DO
	END DO
	
      call fft_xz_to_physical(CS1,S1,0,NY+1)
	
	DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
	      F3(I,K,J) = F3(I,K,J) - ((U2F1(I,K,J+1)+U2F1(I,K,J))*0.5)*
     &			                          S1(I,K,J)
	    END DO
	  END DO
      END DO

C U3: -2*U3F*dU3dX3

	DO J=2,NYM
	  DO K=0,TNKZ
	    DO I=0,NKX
	      CS1(I,K,J) = CIKZ(K)*CU3(I,K,J)
	    END DO
	  END DO
	END DO
	
      call fft_xz_to_physical(CS1,S1,0,NY+1)
	
	DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
	      F3(I,K,J) = F3(I,K,J) - 2.*U3F1(I,K,J)*S1(I,K,J)
	    END DO
	  END DO
	END DO
	
C U3: -U2F*(dU3dX2 + dU2dX3) ??? Interpolations of U2 & U3

	DO J=2,NYM
	  DO K=0,TNKZ
	    DO I=0,NKX
	      CS1(I,K,J) = ( ((0.5*(CU3(I,K,J+1)+CU3(I,K,J)) 
     &			     -0.5*(CU3(I,K,J)+CU3(I,K,J-1)))/DYF(J))
     &			  + (CIKZ(K)*(CU2(I,K,J+1)+CU2(I,K,J))*0.5)  )
	    END DO
	  END DO
	END DO
	
      call fft_xz_to_physical(CS1,S1,0,NY+1)
	
	DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
	      F3(I,K,J) = F3(I,K,J) - ((U2F1(I,K,J+1)+U2F1(I,K,J))/2.0)*
     &					                          S1(I,K,J)
	    END DO
	  END DO
	END DO

C U2: -U1F*(dU2dX1 + dU1dX2) ??? Interpolations of U1

	DO J=2,NYM
	  DO K=0,TNKZ
	    DO I=0,NKX
	      CS1(I,K,J) = (CIKX(I)*CU2(I,K,J) 
     &		     + (CU1(I,K,J)-CU1(I,K,J-1))/DY(J))
	    END DO
	  END DO
	END DO

      call fft_xz_to_physical(CS1,S1,0,NY+1)

	DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
	      F2(I,K,J) = -((U1F1(I,K,J)+U1F1(I,K,J-1))*0.5)*S1(I,K,J)
	    END DO
	  END DO
	END DO

C U2: -2*U2F*dU2dX2		 ??? Interpolations of U2

	DO J=2,NYM
	  DO K=0,TNKZ
	    DO I=0,NKX
	      CS1(I,K,J) = (CU2(I,K,J+1)+CU2(I,K,J))/2.0
	    END DO
	  END DO
	END DO

      call fft_xz_to_physical(CS1,S1,0,NY+1)

	DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
	      F2(I,K,J) = F2(I,K,J) - 2.0*U2F1(I,K,J)*
     $				     ((S1(I,K,J)-S1(I,K,J-1))/DY(J))
	    END DO
	  END DO
	END DO

C U2: -U3F*(dU2dX3 + dU3dX2) ??? Interpolations of U2 & U3

	DO J=2,NYM
	  DO K=0,TNKZ
	    DO I=0,NKX
	      CS1(I,K,J) = (CIKZ(K)*CU2(I,K,J) 
     &		     + (CU3(I,K,J)-CU3(I,K,J-1))/DY(J))
	    END DO
	  END DO
	END DO

      call fft_xz_to_physical(CS1,S1,0,NY+1)

	DO J=2,NYM
        DO K=0,NZM
          DO I=0,NXM
	      F2(I,K,J)=F2(I,K,J)
     %-((U3F1(I,K,J)+U3F1(I,K,J-1))*0.5)*S1(I,K,J)
	    END DO
	  END DO
	END DO


C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!! (NON-LINEAR) FORCING TERMS  !!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C ADD THESE TERMS BEFORE CONVERTING BACK INTO FOURIER SPACE.

C     TO BE DONE

C Convert to Fourier Space
      call FFT_XZ_TO_FOURIER(F1,CF1,0,NY+1)
      call FFT_XZ_TO_FOURIER(F2,CF2,0,NY+1)
      call FFT_XZ_TO_FOURIER(F3,CF3,0,NY+1)
	
C ADD NON-LINEAR TERMS TO THE RHS:
	DO J=2,NYM
	  DO K=0,TNKZ
	    DO I=0,NKX
		CR1(I,K,J)=CR1(I,K,J)+TEMP5*CF1(I,K,J)
		CR3(I,K,J)=CR3(I,K,J)+TEMP5*CF3(I,K,J)
	    END DO
	  END DO
	END DO
	DO J=2,NY
	  DO K=0,TNKZ
	    DO I=0,NKX
		CR2(I,K,J)=CR2(I,K,J)+TEMP5*CF2(I,K,J)
	    END DO
	  END DO
	END DO


C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!! IMPLICIT SOLVES FOR LINEAR TERMS  !!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C The implicit solve in the adjoint case is identical to the state case.
C The majority of the code that follows has been taken directly from
C John Taylor's channel flow solver from Spring 2005.  As updates are
C performed this code will likely be phased out, and replaced by more
C optimized code for this application.


C !!!!!!!!!!!!!!!!!!!!!!!!!!! CU2(I,K,J) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Initialize the matrix to zeros to be used for implicit solves
C Note that the system size is NY+1, but only 1..NY are used
      DO J=0,NY+1
        DO I=0,NKX
          MATL_C(I,J)=0.
          MATD_C(I,J)=0.
          MATU_C(I,J)=0.
        END DO
      END DO

C Build implicit matrix for U2
      DO K=0,TNKZ
        DO J=2,NY
          DO I=0,NKX
            MATL_C(I,J)= -TEMP1/(DYF(J-1)*DY(J))
            MATD_C(I,J)=1.
! Vertical Viscous term
     &         +TEMP1/(DYF(J)*DY(J)) + TEMP1/(DYF(J-1)*DY(J)) 
! Horizontal Viscous terms:
     &         +TEMP1 * (KX2(I)+KZ2(K))
            MATU_C(I,J)= -TEMP1/(DYF(J)*DY(J))
            VEC_C(I,J)=CR2(I,K,J)
          END DO
        END DO

C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Boundary conditions for the adjoint field are non-trivial.  These
C subroutines will need to be worked out to ensure that the are proper. 
!	  CALL APPLY_AD_BC_2_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
!	  CALL APPLY_AD_BC_2_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Tridiagonal Solve for U2(i,:,k)
	  CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY+1,NKX)
        DO J=1,NY+1
          DO I=0,NKX
            CU2(I,K,J)=VEC_C(I,J)
          END DO
        END DO    
      END DO 

C !!!!!!!!!!!!!!!!!!!!!!!!!!! CU1(I,K,J) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Note, here the matrix will be indexed from 1...NY+1 corresponding to U1(0:NY)
      DO J=0,NY+1
        DO I=0,NKX
          MATD_C(I,J)=0.
          MATL_C(I,J)=0.
          MATU_C(I,J)=0.
        END DO
      END DO
      
C Now, solve the implicit equation for U1
      DO K=0,TNKZ
        DO J=2,NYM
          DO I=0,NKX
            MATL_C(I,J)=-TEMP1/(DY(J)*DYF(J))
            MATD_C(I,J)=1.
! Vertical Viscous terms:
     &         +TEMP1*(1./(DY(J+1)*DYF(J))+1./(DY(J)*DYF(J)))
! Horizontal Viscous term:
     &         +TEMP1*(KX2(I)+KZ2(K))
            MATU_C(I,J)=-TEMP1/(DY(J+1)*DYF(J))
            VEC_C(I,J)=CR1(I,K,J)
          END DO
        END DO

C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Boundary conditions for the adjoint field are non-trivial.  These
C subroutines will need to be worked out to ensure that the are proper. 
!	  CALL APPLY_AD_BC_1_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
!	  CALL APPLY_AD_BC_1_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Tridiagonal Solve for U2(:,k,:)
	  CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY+1,NKX)

        DO J=0,NY+1
          DO I=0,NKX
            CU1(I,K,J)=VEC_C(I,J)
          END DO
        END DO
      END DO

C !!!!!!!!!!!!!!!!!!!!!!!!!!! CU3(I,K,J) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Now, solve the implicit equation for U3
      DO K=0,TNKZ
        DO J=2,NYM
          DO I=0,NKX
            MATL_C(I,J)=-TEMP1/(DY(J)*DYF(J))
            MATD_C(I,J)=1.
! Vertical Viscous terms:
     &         +TEMP1*(1./(DY(J+1)*DYF(J))+1./(DY(J)*DYF(J)))
! Horizontal Viscous term:
     &         +TEMP1*(KX2(I)+KZ2(K))
            MATU_C(I,J)=-TEMP1/(DY(J+1)*DYF(J))
            VEC_C(I,J)=CR3(I,K,J)
          END DO
        END DO

C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Boundary conditions for the adjoint field are non-trivial.  These
C subroutines will need to be worked out to ensure that the are proper. 
!	  CALL APPLY_AD_BC_3_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
!	  CALL APPLY_AD_BC_3_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Tridiagonal Solve for U3(i,:,k)
	  CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY+1,NKX)

        DO J=0,NY+1
          DO I=0,NKX
            CU3(I,K,J)=VEC_C(I,J)
          END DO
        END DO
      END DO
C -- Done getting CU1, CU2, CU3 at new RK Step --


C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!! DIVERGENCE REMOVAL  !!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Similar to the implicit solves.  This code is very similar to the code
C written by John Taylor in Spring 2005 for the channel flow case.  The
C subroutines below are based on equivalent version that John applies 
C in his code.

C Begin second step of the Fractional Step algorithm, making u divergence free
C The following subroutine projects Uhat onto divergence free space

      CALL REM_DIV_ADJ_CHAN

C Now, phi is stored in CR1, use this to update the pressure field
C Note, here we divide by H_BAR since it was absorbed into PHI in REM_DIV
      DO J=1,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CP(I,K,J)=CP(I,K,J)+CR1(I,K,J)/TEMP4
          END DO
        END DO
      END DO

      RETURN
      END
	
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_ADJ_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

C This code solves for the appropriate update which makes this flow
C divergence free.  In effect the adjoint state is being projected onto
C a divergence free manifold.  This code is very similar to John 
C Taylor's REM_DIV_CHAN code.  A few minor changes have been made to 
C adapt the solve for the adjoint case.

C Compute varphi, store in variable CR1.
C Solves for phi in computational space
C H_BAR has been absorbed into PHI, so we are solving for H_BAR*PHI

	INCLUDE 'header_batch'
	INTEGER I,J,K
 
C First, Initialize the matrix components
      DO J=0,NY+1
        DO I=0,NKX
          MATL_C(I,J)=0.
          MATD_C(I,J)=1.
          MATU_C(I,J)=0.
          VEC_C(I,J)=(0.,0.)
        END DO
      END DO

C The 2d FFT of Ui should have been taken and stored in CUi
C Solving for phi amounts to solving a tridiagonal system
C First, construct the system to be solved
      DO K=0,TNKZ
        DO J=1,NY
          DO I=0,NKX
            MATL_C(I,J)=1./(DY(J)*DYF(J))
            MATD_C(I,J)=-KX2(I)-KZ2(K)
     &         -1./(DY(J+1)*DYF(J))-1./(DY(J)*DYF(J))
            MATU_C(I,J)=1./(DY(J+1)*DYF(J))
          END DO
        END DO

C Now, create the RHS vector
        DO J=1,NY         
          DO I=0,NKX
            VEC_C(I,J)=(CIKX(I)*CU1(I,K,J) 
     &               + (CU2(I,K,J+1)-CU2(I,K,J))/DYF(J) 
     &               +  CIKZ(K)*CU3(I,K,J))
	    END DO
	  END DO

        DO I=0,NKX
          IF ((K.EQ.0).AND.(I.EQ.0)) THEN
C Use homogeneous dirichlet BCS for kx=kz=0 component at bottom wall
C Otherwise the matrix will be singular
            MATL_C(I,1)=0. 
            MATD_C(I,1)=1.
            MATU_C(I,1)=0.
            VEC_C(I,1)=(0.,0.)

            MATL_C(I,NY)=1.
            MATD_C(I,NY)=-1.
            MATU_C(I,NY)=0.
            VEC_C(I,NY)=(0.,0.)
          ELSE
C Use Dirichlet boundary conditions, dp/dz=0 at walls
            MATL_C(I,1)=0.
            MATD_C(I,1)=1.
            MATU_C(I,1)=-1.
            VEC_C(I,1)=(0.,0.)

            MATL_C(I,NY)=1.
            MATD_C(I,NY)=-1.
            MATU_C(I,NY)=0.
            VEC_C(I,NY)=(0.,0.)
          END IF
        END DO

C Now solve the tridiagonal system for phi, store in CR1
        CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY,NKX)

C REPOPULATE CR1 WITH THE UPDATE:
        DO J=1,NY
          DO I=0,NKX
            CR1(I,K,J)=VEC_C(I,J)
          END DO
        END DO

      END DO

C Now, Solve for CUi, the divergenceless velocity field
      DO J=1,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CU1(I,K,J)=CU1(I,K,J)-CIKX(I)*CR1(I,K,J)
            CU3(I,K,J)=CU3(I,K,J)-CIKZ(K)*CR1(I,K,J)           
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CU2(I,K,J)=CU2(I,K,J)-(CR1(I,K,J)
     &                                -CR1(I,K,J-1))/DY(J)
          END DO
        END DO
      END DO

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_AD_BC_1_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header_batch'
      INTEGER I

C Bottom Wall:
      IF (U_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,0)=0. 
          MATD_C(I,0)=1.
          MATU_C(I,0)=0.                   
          VEC_C(I,0)=0.

          MATL_C(I,1)=0. 
          MATD_C(I,1)=1.
          MATU_C(I,1)=0.
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
          VEC_C(I,1)=U_BC_YMIN_C1 
        END DO
      ELSE
C Neumann
        DO I=0,NKX
          MATL_C(I,0)=0.
          MATD_C(I,0)=-1.
          MATU_C(I,0)=1.
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
          VEC_C(I,0)=DY(1)*U_BC_YMIN_C1
        END DO
      END IF

      RETURN 
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_AD_BC_1_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header_batch'
      INTEGER I

C Top wall
      IF (U_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          VEC_C(I,NY+1)=0.

          MATL_C(I,NY)=0.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
          VEC_C(I,NY)=U_BC_YMAX_C1
        END DO
      ELSE
C Neumann
        DO I=0,NKX
          MATL_C(I,NY+1)=-1.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
          VEC_C(I,NY+1)=DY(NY+1)*U_BC_YMAX_C1
        END DO      
      END IF

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_AD_BC_2_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header_batch'
      INTEGER I

C Bottom Wall:
      IF (V_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,1)=0.d0 
          MATD_C(I,1)=1.d0
          MATU_C(I,1)=0.d0                   
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
          VEC_C(I,1)=V_BC_YMIN_C1 

          MATL_C(I,2)=0.d0 
          MATD_C(I,2)=1.d0
          MATU_C(I,2)=0.d0                   
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
          VEC_C(I,2)=V_BC_YMIN_C1 
        END DO
      ELSE IF (V_BC_YMIN.EQ.1) THEN
C Neumann
        DO I=0,NKX
          MATD_C(I,1)=-1.d0
          MATU_C(I,1)=1.d0
          MATL_C(I,1)=0.d0
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
         VEC_C(I,1)=DYF(1)*V_BC_YMIN_C1
        END DO
      END IF

C The following is only a placeholder, this row is used for U1 and U3
      DO I=0,NKX
        MATL_C(I,0) = 0.
        MATD_C(I,0) = 1.
        MATU_C(I,0) = 0.
        VEC_C(I,0) = 0.
      END DO

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_AD_BC_2_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header_batch'
      INTEGER I
C Top wall
      IF (V_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
          VEC_C(I,NY+1)=V_BC_YMAX_C1
          
          MATL_C(I,NY)=0.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
          VEC_C(I,NY)=V_BC_YMAX_C1
        END DO
      ELSE IF (V_BC_YMAX.EQ.1) THEN
C Neumann
        DO I=0,NKX
          MATL_C(I,NY+1)=-1.
          MATD_C(I,NY+1)=1.
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
          VEC_C(I,NY+1)=DYF(NY)*V_BC_YMAX_C1
        END DO      
      END IF
	
      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_AD_BC_3_LOWER_C(MATL_C,MATD_C,MATU_C,VEC_C)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header_batch'
      INTEGER I

C Bottom Wall:
      IF (W_BC_YMIN.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,0)=0. 
          MATD_C(I,0)=1.
          MATU_C(I,0)=0.                   
          VEC_C(I,0)=0.

          MATL_C(I,1)=0. 
          MATD_C(I,1)=1.
          MATU_C(I,1)=0.                   
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
          VEC_C(I,1)=W_BC_YMIN_C1
        END DO
      ELSE
C Neumann
        DO I=0,NKX
          MATL_C(I,0)=0.
          MATD_C(I,0)=-1.
          MATU_C(I,0)=1.
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
          VEC_C(I,0)=DY(1)*W_BC_YMIN_C1
        END DO
      END IF

      RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_AD_BC_3_UPPER_C(MATL_C,MATD_C,MATU_C,VEC_C)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header_batch'
      INTEGER I

C Top wall
      IF (W_BC_YMAX.EQ.0) THEN
C Dirichlet
        DO I=0,NKX
          MATL_C(I,NY+1)=0.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
          VEC_C(I,NY+1)=0.

          MATL_C(I,NY)=0.
          MATD_C(I,NY)=1.
          MATU_C(I,NY)=0.
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
          VEC_C(I,NY)=W_BC_YMAX_C1
        END DO
      ELSE
C Neumann
        DO I=0,NKX
          MATL_C(I,NY+1)=-1.
          MATD_C(I,NY+1)=1.
          MATU_C(I,NY+1)=0.
C CHANGE THE VALUE BELOW FOR ADJOINT BC'S
          VEC_C(I,NY+1)=DY(NY+1)*W_BC_YMAX_C1
        END DO      
      END IF

      RETURN
      END
	
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE output_adjoint(plane,vloc,noise)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This code extracts a plane of measurements for the sensors located at
! positions defined by vloc.  If these measurements are being extracted
! from a simulation then noise can be added to the values by choosing
! noise=1.  If these measurements are being extracted from the estimate,
! then no noise should be added, and noise should be inputed as 0 (zero)
	include 'header_batch'
! Local Variables for output
      INTEGER m, v, i, k, noise
	REAL*8 zbqlnor, mu, vloc(1:3,1:nm), plane(1:nv,1:nm)
! Measurement variables passed out by find_meas_loc:
	REAL*8  abs_err
	INTEGER loc
	
!	Statistics of Noise
	mu = 0.


	if (num_per_dir.eq.2) then
	
	plane = 0.0
	
	do m=1,NM
	  select case (meas_type_adjoint(m))

C	    Pressure Measurement Extraction
	    case ('P')
	      if (verbosity.gt.2) then
	      write(6,*) 'Getting Pressure'
		end if
	      do v=1,NV
		  if (verbosity.gt.2) then
		    write(6,*) 'Sensor Location: ', vloc(2,v)
		  end if
		  call find_meas_loc_y(vloc(2,v),loc,abs_err,1)
		  do i=0,nkx
		    do k=0,tnkz
		      plane(v,m) = plane(v,m) +
     &		 ((abs_err*CP(i,k,loc+1)+(1.0-abs_err)*CP(i,k,loc)))*
     &         (exp(-(0.5*(ave_var**2.)*(kx2(i)+kz2(k)))
     &        -(cikx(i)*(LX-vloc(1,v)) + cikz(k)*(LZ-vloc(3,v)))))
     			if (i.gt.0) then
		  plane(v,m) = plane(v,m) + 
     &(((abs_err*conjg(CP(i,k,loc+1)))+(1.-abs_err)*conjg(CP(i,k,loc)))*
     &        (exp(-(0.5*(ave_var**2.)*(kx2(i)+kz2(k)))
     &       +(cikx(i)*(LX-vloc(1,v)) + cikz(k)*(LZ-vloc(3,v))))))
			end if
		    end do
		  end do
		end do
		
C	    Shear (x-direction) Measurement Extraction
	    case('SH1')
	      if (verbosity.gt.2) then
	        write(6,*) 'Getting Wall-Shear (X-Direction)'
		end if
	      do v=1,NV
		  call find_meas_loc_y(vloc(2,v),loc,abs_err,1)
		  do i=0,nkx
		    do k=0,tnkz
		      plane(v,m) = plane(v,m) +
     &		((abs_err*CU1(i,k,loc+1)+(1.0-abs_err)*CU1(i,k,loc)))*
     &         (exp(-(0.5*(ave_var**2.)*(kx2(i)+kz2(k)))
     &              -(cikx(i)*(LX-vloc(1,v)) + cikz(k)*(LZ-vloc(3,v)))))
     			if (i.gt.0) then
		  plane(v,m) = plane(v,m) +(((abs_err*conjg(CU1(i,k,loc+1)))
     &                              + (1.-abs_err)*conjg(CU1(i,k,loc)))*
     &        (exp(-(0.5*(ave_var**2.)*(kx2(i)+kz2(k)))
     &             +(cikx(i)*(LX-vloc(1,v)) + cikz(k)*(LZ-vloc(3,v))))))
			end if
		    end do
		  end do
		  plane(v,m) = nu*plane(v,m)/dy(2)
		end do
	    
C	    Shear (z-direction) Measurement Extraction
	    case('SH3')
	      if (verbosity.gt.2) then
	        write(6,*) 'Getting Wall-Shear (Z-Direction)'
	      end if
		do v=1,NV
		  call find_meas_loc_y(vloc(2,v),loc,abs_err,1)
		  do i=0,nkx
		    do k=0,tnkz
		      plane(v,m) = plane(v,m) +
     &		((abs_err*CU3(i,k,loc+1)+(1.0-abs_err)*CU3(i,k,loc)))*
     &         (exp(-(0.5*(ave_var**2.)*(kx2(i)+kz2(k)))
     &              -(cikx(i)*(LX-vloc(1,v)) + cikz(k)*(LZ-vloc(3,v)))))
     			if (i.gt.0) then
		  plane(v,m) = plane(v,m)+ (((abs_err*conjg(CU3(i,k,loc+1)))
     &                    + (1.-abs_err)*conjg(CU3(i,k,loc)))*
     &        (exp(-(0.5*(ave_var**2.)*(kx2(i)+kz2(k)))
     &             +(cikx(i)*(LX-vloc(1,v)) + cikz(k)*(LZ-vloc(3,v))))))
			end if
		    end do
		  end do
		  plane(v,m) = nu*plane(v,m)/dy(2)
		end do
	    case default
	    write(6,*) 'You Suck: Measurement Type not supported.'
	  end select
	
C	  If Noise is to be added.  DO IT HERE!!!
	  if (noise.eq.1) then
	    do v=1,nv
		plane(v,m) = plane(v,m) + zbqlnor(mu,sqrt(meas_var(m)))
	    end do
	  end if
	  
	end do
	  
      else
	  write(6,*) 'You suck: Not using channel case (output_adjoint)'
      end if
	
      RETURN
      END
	
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE find_meas_loc_y(position,loc,abs_err,kind)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine maps the y-position of a vehicle in vpos to the
! j-plane (plane-plane) directly below the point.  It outputs the plane
! integer number and a value, bracketed by [0,1), describing where the
! point is between that j-plane and the plane directly above it.  Also,
! two different finds are needed: one for the base grid, and one for the
! fractional grid
	include 'header_batch'
      
C	Local Variables:
	INTEGER kind
	real*8 error, temp

! Measurement variables passed out by find_meas_loc:
	REAL*8  position, abs_err
	INTEGER loc

	temp = LY
	loc = ny+1
	
C	There are two types of measurement location.  Type ZERO corresponds
C	values that are defined on the base grid.  Type ONE corresponds to
C	values that are defined on the fractional grid.
	select case (kind)
	
C	  ! Base Grid
	  case (0)
	    write(6,*)
     &       'You suck: Base Grid Code not defined (find_meas_loc_y).'
	  
C	  ! Fractional Grid
	  case (1)
	    do while (gyf(loc).gt.position)
		loc = loc-1
c	      write(6,*) 'GY: ', gyf(loc),'POS: ', pos, ' Loc: ', loc
	    end do
	    abs_err = (position-gyf(loc))/(gyf(loc+1)-gyf(loc))
	    if (verbosity.gt.2) then
	      write(6,*) 'Lower Plane: ', loc,
     &	 'Lower Plane Height: ', gyf(loc),
     &	 'Position: ', position,
     &	 'Upper Plane Height: ', gyf(loc+1),
     &	 'abs_err: ', abs_err,
     &       'Ratio Check: ', (position-gyf(loc))/(gyf(loc+1)-gyf(loc)),
     &       'Grid Width', dy(loc+1)
	    end if
	    
C	  ! Otherwise
	  case default
	    write(6,*) 'You suck: Define grid type (find_meas_loc_y).'
	
	end select



      RETURN
      END
	
	
	
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine test_brent
 	real*8 ax,bx,cx,fa,fb,fc,freq,xmin

	  
c THE CODE BELOW TESTS TO MAKE SURE THAT MNBRAK AND BREANT ARE WORKING !	
	freq=1000;
	pi = 4.0*atan(1.0)
	ax = 10
	bx = 100
	call cost_function(fa,ax)
	call cost_function(fb,bx)
	write(*,*) ax*(pi/freq), fa
	write(*,*) bx*(pi/freq), fb
	call mnbrak(ax,bx,cx,fa,fb,fc)
	write(*,*) ax*(pi/freq), fa
	write(*,*) bx*(pi/freq), fb
	write(*,*) cx*(pi/freq), fc
	call brent(ax,bx,cx,1.e-10,xmin)
	write(*,*) xmin*(pi/freq)
	  
	return
	end
	
	
	
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine test_measurements
	  
c !!!!!!!!!!!!!!!! TEST CODE TO CHECK MEASUREMENTS !!!!!!!!!!!!!!!!!!!!
c To check that the measurements are extracted properly use the following 
c lines.  Make sure to remove the scaler constants for shear 
c in output_adjoint!
c	  if (mod(TIME_STEP,100).eq.0) then
c	  call output_adjoint(y(1:nv,1:nm,1),veh_pos(1:3,1:nv,1),0)
c        call fft_xz_to_physical(CP,P,0,NY+1)
c        call fft_xz_to_physical(CU1,U1,0,NY+1)
c        call fft_xz_to_physical(CU3,U3,0,NY+1)
c	    write(6,*) 'Actaul Measurements: ', P(15,4,1), P(19,16,1),
c     &	U1(15,4,1), U1(19,16,1), U3(15,4,1), U3(19,16,1)
c	  v = 1
c	  do i=0,nxm
c	    do j=0,nzm
c		write(6,*) 'Error: ', v, y(v,1,1) - P(i,j,2)
c     &				     , y(v,2,1) - U1(i,j,2)
c     &				     , y(v,3,1) - U3(i,j,2)
c		v = v+1
c	    end do
c	  end do
c	  do i=0,nxm
c	    do j=0,nzm
c		write(6,*) 'Error: ', v, y(v,1,1) - P(i,j,nym)
c     &				     , y(v,2,1) - U1(i,j,nym)
c     &				     , y(v,3,1) - U3(i,j,nym)
c		v = v+1
c	    end do
c	  end do
c        call fft_xz_to_fourier(P,CP,0,NY+1)
c        call fft_xz_to_fourier(U1,CU1,0,NY+1)
c        call fft_xz_to_fourier(U3,CU3,0,NY+1)
c	  end if
c !!!!!!!!!!!!!!!!!! END TEST CODE FOR MEASUREMENTS !!!!!!!!!!!!!!!!!!!!
	  
	return
	end	
	
	
	
	
	
	