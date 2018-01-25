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
      REAL*8 TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6

C Define the constants that are used in the time-stepping
C For reference, see Numerical Renaissance
      TEMP1=NU * H_BAR(RK_STEP) / 2.0
      TEMP2=BETA_BAR(RK_STEP)*H_BAR(RK_STEP)
      TEMP3=ZETA_BAR(RK_STEP) * H_BAR(RK_STEP)
      TEMP4=H_BAR(RK_STEP)
	TEMP5=0.0

	if ((FLAVOR.eq.'Batch') .AND. (delta_t.lt.0.0)) then
	  call backwards_constants(TEMP1,TEMP2,TEMP3,TEMP4,TEMP5)
	end if

C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!  INCLUDE PREVIOUS ADJOINT STATES  !!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
C Start with the explicit part of the Crank-Nicolson viscous term and
C  the pressure gradient treated with Explicit Euler:
            TEMP5=1+TEMP1*(KX2(I)+KY2(J)+KZ2(K))
		
C The following hook should be used for regularizing the NS equations.
C If marches are backwards-in-time (ie delta_t<0) then it will be used.
C		IF (DELTA_T.LT.0) CALL quasi_rev_per(TEMP5,i,j,k)
		
            CR1(I,K,J)=TEMP5*CU1S(I,K,J)-TEMP4*(CIKX(I)*CPS(I,K,J))
            CR2(I,K,J)=TEMP5*CU2S(I,K,J)-TEMP4*(CIKY(J)*CPS(I,K,J))
            CR3(I,K,J)=TEMP5*CU3S(I,K,J)-TEMP4*(CIKZ(K)*CPS(I,K,J))
C For each scalar, start with the explict part of the Crank-Nicolson
C diffusive term for each scalar
            DO N=1,N_TH
              TEMP6=1+(TEMP1/PR(N))*(KX2(I)+KY2(J)+KZ2(K))
     &               -REACTION(N)*TEMP1
              CRTH(I,K,J,N)=TEMP6*CTHS(I,K,J,N)
            END DO
          END DO
        END DO
        IF (RK_STEP .GT. 1) THEN
          DO K=0,TNKZ
            DO I=0,NKX
C Add the term: ZETA_BAR(RK_STEP)*R(U(RK_STEP-1))
              CR1(I,K,J)=CR1(I,K,J)+TEMP3*CF1S(I,K,J)
              CR2(I,K,J)=CR2(I,K,J)+TEMP3*CF2S(I,K,J)
              CR3(I,K,J)=CR3(I,K,J)+TEMP3*CF3S(I,K,J)
C Do the same for each scalar:
              DO N=1,N_TH
                CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP3*CFTHS(I,K,J,N)
              END DO
            END DO
          END DO
        END IF
      END DO


C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!  NON-LINEAR TERMS  !!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Convert flow velocities to physical space, because non-linear products
C in the adjoint equation need to be calculated in physical space.
      CALL FFT_XZY_TO_PHYSICAL(CU1,U1)
      CALL FFT_XZY_TO_PHYSICAL(CU2,U2)
      CALL FFT_XZY_TO_PHYSICAL(CU3,U3)
      DO N=1,N_TH 
        CALL FFT_XZY_TO_PHYSICAL(CTHS(0,0,0,N),THS(0,0,0,N))
      END DO


C !!!!!!!!!!!!!!!!!!!!!  PASSIVE SCALAR TERMS  !!!!!!!!!!!!!!!!!!!!!!!!!
C Compute the nonlinear terms for the passive scalar equation
C Do this before the nonlinear momentum terms to use Fi as a working
C  array before using it for the momentum equation.
c      IF (TRUCK) CALL INIT_FORCING

      DO N=1,N_TH
        DO J=0,NYM
          DO K=0,NZM
            DO I=0,NXM
              F1S(I,K,J)=U1(I,K,J)*THS(I,K,J,N)
              F2S(I,K,J)=U2(I,K,J)*THS(I,K,J,N)
              F3S(I,K,J)=U3(I,K,J)*THS(I,K,J,N)
            END DO
          END DO
        END DO
        CALL FFT_XZY_TO_FOURIER(F1S,CF1S)
        CALL FFT_XZY_TO_FOURIER(F2S,CF2S)
        CALL FFT_XZY_TO_FOURIER(F3S,CF3S)
        DO J=0,TNKY
          DO K=0,TNKZ
            DO I=0,NKX
              CFTHS(I,K,J,N)=CFTHS(I,K,J,N)-CIKX(I)*CF1S(I,K,J)
     &                      -CIKY(J)*CF2S(I,K,J)
     &                      -CIKZ(K)*CF3S(I,K,J)
     
C Add the forcing due to the moving source
c              IF (TRUCK) CALL FORCING(I,K,J)
     
C Add R-K terms for the TH equation to the RHS
              CRTH(I,K,J,N)=CRTH(I,K,J,N)+TEMP2*CFTHS(I,K,J,N)
            END DO
          END DO
        END DO
      END DO
C The RHS vector for the TH equation is now ready


C !!!!!!!!!!!!!!!!!!!!  MOMENTUM EQUATION TERMS  !!!!!!!!!!!!!!!!!!!!!!!
C Compute the nonlinear terms for the momentum equations

      DO j=0,TNKY
        DO k=0,TNKZ
          DO i=0,NKX
            CF1S(i,k,j)= CIKX(i)*CU1S(i,k,j)
            CF2S(i,k,j)= CIKY(j)*CU1S(i,k,j) + CIKX(i)*CU2S(i,k,j)
            CF3S(i,k,j)= CIKZ(k)*CU1S(i,k,j) + CIKX(i)*CU3S(i,k,j)
		CS1(i,k,j)= CIKY(j)*CU2S(i,k,j)
          END DO
        END DO
       END DO
      CALL FFT_XZY_TO_PHYSICAL(CF1S,F1S)
      CALL FFT_XZY_TO_PHYSICAL(CF2S,F2S)
      CALL FFT_XZY_TO_PHYSICAL(CF3S,F3S)
      CALL FFT_XZY_TO_PHYSICAL(CS1,S1)
C Here we start constructing the R-K terms in CFi
C Note, that the order of the following operations are important

      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            F1S(i,k,j) = - 2.0*F1S(i,k,j)*U1(i,k,j)
     &                   -     F2S(i,k,j)*U2(i,k,j)
     &                   -     F3S(i,k,j)*U3(i,k,j)
            F2S(i,k,j) = -     F2S(i,k,j)*U1(I,K,J)
     &                   -  2.0*S1(i,k,j)*U2(i,k,j)
            F3S(i,k,j) = -     F3S(i,k,j)*U1(i,k,j)
          END DO
        END DO
      END DO


      DO j=0,TNKY
        DO k=0,TNKZ
          DO i=0,NKX
		CS1(i,k,j)= CIKY(j)*CU3S(i,k,j) + CIKZ(k)*CU2S(i,k,j)
          END DO
        END DO
       END DO
      CALL FFT_XZY_TO_PHYSICAL(CS1,S1)

      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            F2S(i,k,j) = F2S(i,k,j) - S1(i,k,j)*U3(I,K,J)
            F3S(i,k,j) = F3S(i,k,j) - S1(i,k,j)*U2(i,k,j)
          END DO
        END DO
      END DO

      DO j=0,TNKY
        DO k=0,TNKZ
          DO i=0,NKX
		CS1(i,k,j)= CIKZ(k)*CU3S(i,k,j)
          END DO
        END DO
       END DO
      CALL FFT_XZY_TO_PHYSICAL(CS1,S1)

      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            F3S(i,k,j) = F3S(i,k,j) - 2.0*S1(i,k,j)*U3(i,k,j)
          END DO
        END DO
      END DO

c	Scaler Forcing: The scaler naturally forces the adjoint equations
c	for the velocities in each direction.
	DO N=1,N_TH

	  DO j=0,TNKY
          DO k=0,TNKZ
            DO i=0,NKX
		  CS1(i,k,j)= CIKX(i)*CTH(i,k,j,N)
            END DO
          END DO
         END DO
        CALL FFT_XZY_TO_PHYSICAL(CS1,S1)
	  
        DO J=0,NYM
          DO K=0,NZM
            DO I=0,NXM
              F1S(i,k,j) = F1S(i,k,j) + S1(i,k,j)*THS(i,k,j,N)
            END DO
          END DO
        END DO
	  
        DO j=0,TNKY
          DO k=0,TNKZ
            DO i=0,NKX
		  CS1(i,k,j)= CIKY(j)*CTH(i,k,j,N)
            END DO
          END DO
         END DO
        CALL FFT_XZY_TO_PHYSICAL(CS1,S1)
	  
        DO J=0,NYM
          DO K=0,NZM
            DO I=0,NXM
              F2S(i,k,j) = F2S(i,k,j) + S1(i,k,j)*THS(i,k,j,N)
            END DO
          END DO
        END DO	  
	  
        DO j=0,TNKY
          DO k=0,TNKZ
            DO i=0,NKX
		  CS1(i,k,j)= CIKZ(k)*CTH(i,k,j,N)
            END DO
          END DO
         END DO
        CALL FFT_XZY_TO_PHYSICAL(CS1,S1)
	  
        DO J=0,NYM
          DO K=0,NZM
            DO I=0,NXM
              F3S(i,k,j) = F3S(i,k,j) + S1(i,k,j)*THS(i,k,j,N)
            END DO
          END DO
        END DO

      END DO


C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!! (NON-LINEAR) FORCING TERMS  !!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Rather than forcing the RHS of this equation, we are performing discrete
c updates to adjoint at the same time as we remove directions in the backwards
c direction.  This will be done in the main body of the code.

C Convert to Fourier Space
      call FFT_XZY_TO_FOURIER(F1S,CF1S)
      call FFT_XZY_TO_FOURIER(F2S,CF2S)
      call FFT_XZY_TO_FOURIER(F3S,CF3S)
	
C ADD NON-LINEAR TERMS TO THE RHS:
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CR1(I,K,J)=CR1(I,K,J)+TEMP2*CF1(I,K,J)
            CR2(I,K,J)=CR2(I,K,J)+TEMP2*CF2(I,K,J)
            CR3(I,K,J)=CR3(I,K,J)+TEMP2*CF3(I,K,J)
          END DO
        END DO
      END DO
C Computation of CRi complete.


C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!! IMPLICIT SOLVES FOR LINEAR TERMS  !!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Now solve the implicit system for the intermediate field.
C (In the fully-periodic case, this is easy!)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            TEMP5=1-TEMP1*(KX2(I)+KY2(J)+KZ2(K))

C The following hook should be used for regularizing the NS equations.
C If marches are backwards-in-time (ie delta_t<0) then it will be used.		

		CU1(I,K,J)=CR1(I,K,J)/TEMP5
            CU2(I,K,J)=CR2(I,K,J)/TEMP5
            CU3(I,K,J)=CR3(I,K,J)/TEMP5
            DO N=1,N_TH
              TEMP6=1-(TEMP1/PR(N))*(KX2(I)+KY2(J)+KZ2(K)) 
     &		   + REACTION(N)*TEMP1
              CTHS(I,K,J,N)=CRTH(I,K,J,N)/TEMP6
            END DO
          END DO
        END DO
      END DO
C First step of the Fractional Step algorithm complete.


C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!! DIVERGENCE REMOVAL  !!!!!!!!!!!!!!!!!!!!!!!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C Begin second step of the Fractional Step algorithm, making u divergence free
C The following subroutine projects Uhat onto divergence free space

      CALL REM_DIV_ADJ_PER

C Now, phi is stored in CR1, use this to update the pressure field
C Note, here we divide by H_BAR since it was absorbed into PHI in REM_DIV
      DO J=1,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CP(I,K,J)=CP(I,K,J) - CR1(I,K,J)/TEMP4
          END DO
        END DO
      END DO

      RETURN
      END
	
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE REM_DIV_ADJ_PER
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
      REAL*8  TEMP5

C Compute phi, store in the variable CR1.
C Note the coefficient H_BAR is absorbed into phi.
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            TEMP5 = -(KX2(I)+KY2(J)+KZ2(K)+EPS)
            CR1(I,K,J)=(CIKX(I)*CU1S(I,K,J)+CIKY(J)*CU2S(I,K,J)+
     *                  CIKZ(K)*CU3S(I,K,J))/TEMP5
          END DO
        END DO
C Then update the CUi to make velocity field divergence-free.
        DO K=0,TNKZ
          DO I=0,NKX
            CU1S(I,K,J)=CU1S(I,K,J)-CIKX(I)*CR1(I,K,J)
            CU2S(I,K,J)=CU2S(I,K,J)-CIKY(J)*CR1(I,K,J)
            CU3S(I,K,J)=CU3S(I,K,J)-CIKZ(K)*CR1(I,K,J)
          END DO
        END DO
      END DO

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
	 write(6,*) 'Adjoint Output: Not using channel cases.'
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
	
	
	
	
	
	