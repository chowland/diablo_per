C******************************************************************************|
C name.f -> DNS In A Box, Laptop Optimized                          VERSION 0.9
C
C This Fortran 77 code estimates incompressible flow for a box.
C
C Primative variables (u,v,w,p) are used, and continuity is enforced with a
C fractional step algorithm.
C
C SPATIAL DERIVATIVES:
C   0, 1, 2, or 3 directions are taken to be periodic and handled spectrally
C   (these cases are referred to as the "periodic", "channel", "duct", and
C    "cavity" cases respectively).
C   The remaining directions are taken to be bounded by walls and handled with
C   momentum- and energy-conserving second-order central finite differences.
C
C TIME ADVANCEMENT
C   Two main approaches are implemented:
C     1. RKW3 on nonlinear terms and CN on viscous terms over each RK substep.
C     2. RKW3 on y-derivative terms and CN on other terms over each RK substep.
C
C The emphasis in this introductory code is on code simplicity:
C   -> All variables are in core.
C   -> The code is not explicitly designed for use with either MPI or SMP.
C   -> Overindexing is not used.
C A few simple high-performance programming constructs are used:
C   -> The inner 2 loops are broken out in such a way as to enable out-of-order
C      execution of these loops as much as possible, thereby leveraging
C      vector and superscalar CPU architectures.
C   -> The outer loops are fairly long (including as many operations as
C      possible inside on a single J plane of data) in order to make effective
C      use of cache.
C Multiple time advancement algorithms are implemented for the periodic,
C channel, duct, and cavity cases in order to compare their efficiency for
C various flows on different computational architectures.  In a few of the
C algorithms, the code is broken out fully into a PASS1/PASS2 architecture
C to maximize the efficient use of cache.
C
C This code was developed as a joint project between the University of 
C California - San Diego and the National Secutiy Education Center - Los
C Alamos National Laboratories.  The origianl code was developed jointly
C by undergraduate and graduate students at UCSD in MAE 223 (CFD), taught by
C Thomas Bewley, at UC San Diego (spring of 2001, 2005, 2007).
C
C Primary contributions follow:
C Thomas Bewley was the chief software architect
C John R. Taylor wrote the channel flow solvers
C Chris Colburn wrote adjoint solvers and minimization algorithms
C******************************************************************************|
C
C This code is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the
C Free Software Foundation; either version 2 of the License, or (at your
C option) any later version. This code is distributed in the hope that it
C will be useful, but WITHOUT ANY WARRANTY; without even the implied
C warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C GNU General Public License for more details. You should have received a
C copy of the GNU General Public License along with this code; if not,
C write to the Free Software Foundation, Inc., 59 Temple Place - Suite
C 330, Boston, MA 02111-1307, USA.
C
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      PROGRAM TEST_CODES
      WRITE(6,*)
      WRITE(6,*) '             ***************************************'
      WRITE(6,*) '             ******   WELCOME TO THE BATCH    ******'
      WRITE(6,*) '             ******  PROCESSING UNIT OF EnVE  ******'
      WRITE(6,*) '             ***************************************'
      WRITE(6,*)
	
	call REGULARIZATION
c	call BATCH
	
C Below are subroutines written to make troubleshooting code streamlined
C 

c	TEST PENTA_SOLVE
c	call test_penta
C	call test_septa

C If you want to test the minimization tool.
c	call test_brent
	
      WRITE(6,*)
      WRITE(6,*) '        ****** Hello world!  Have a nice day! ******'
      WRITE(6,*)
      
	return
	END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine BATCH
C This code is designed to find statistical averages of L2 norms of the
C error in a backwards march of the Navier Stokes Equations regularized
C with a Mized Time-Space Derivative.  The code DOES NOT FIND GRAMMER
C MISTAKES, like the fragment before this sentance.
	INCLUDE 'header_batch'
      INTEGER ind, INITIALIZE_TIME, LAST_TIME, MAXWIN, i, k, j
	complex*16 store_CU1(0:NX,0:NZ,0:NY+1)
     &         , store_CU2(0:NX,0:NZ,0:NY+1)
     &         , store_CU3(0:NX,0:NZ,0:NY+1)
     	real*8 sum1(0:NX,0:NZ), sum2(0:NX,0:NZ), sum3(0:NX,0:NZ)
	real*8 scaler1, scaler3

C	Change INITIALIZE so that the appropriate estimation field is 
C	being read into name.f
      CALL INITIALIZE
	call init_march

! Initialize START_TIME for run timing
      CALL DATE_AND_TIME (VALUES=TIME_ARRAY)
      START_TIME=TIME_ARRAY(5)*3600+TIME_ARRAY(6)*60
     &   +TIME_ARRAY(7)+TIME_ARRAY(8)*0.001

C A flag to determine if we are considering the first time-step
      FIRST_TIME=.TRUE.

C INITIALIZE FLOW TO GET A SMOOTH TRAJECTORY AND POPULATE MEASUREMENTS
	INITIALIZE_TIME=500
	GUSTS = .TRUE.
	DO TIME_STEP = 1,INITIALIZE_TIME
        WRITE(6,*) 'Now ending TIME_STEP = ',TIME_STEP
C	  Forward Time Evolution of State
        DO RK_STEP=1,3
          IF (NUM_PER_DIR.EQ.3) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_PER_1
            IF (TIME_AD_METH.EQ.2) CALL RK_PER_2            
          ELSEIF (NUM_PER_DIR.EQ.2) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CHAN_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CHAN_2            
          ELSEIF (NUM_PER_DIR.EQ.1) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_DUCT_1
            IF (TIME_AD_METH.EQ.2) CALL RK_DUCT_2            
          ELSEIF (NUM_PER_DIR.EQ.0) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CAV_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CAV_2            
          END IF
        END DO

C	  Extract Measurements when Appropriate
c	  if (mod(TIME_STEP,smpl_freq).eq.0) then
c	    ind = mod(TIME_STEP/smpl_freq,THL)
c	    call poisson_p_chan
c	    call output_adjoint(y(1:nv,1:nm,ind),veh_pos(1:3,1:nv,ind),0)
c	  end if	
	  
	  if (mod(TIME_STEP,100).eq.0) then
	    do j = 0,ny+1
		write(*,'(I4,F12.6,E12.3,E12.3)') j, real(CU1(0,0,j)),
     &                   real(CU2(0,0,j)), real(CU3(0,0,j))
	    end do
	  end if

	END DO
	
	write(6,*) 'Flow has been initialized!'
	GUSTS = .FALSE.
	LAST_TIME = TIME_STEP
	MAXWIN=10

      DO TIME_STEP = TIME_STEP, TIME_STEP+MAXWIN-1

C	  Forward Time Evolution of State
        DO RK_STEP=1,3
          IF (NUM_PER_DIR.EQ.3) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_PER_1
            IF (TIME_AD_METH.EQ.2) CALL RK_PER_2            
          ELSEIF (NUM_PER_DIR.EQ.2) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CHAN_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CHAN_2            
          ELSEIF (NUM_PER_DIR.EQ.1) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_DUCT_1
            IF (TIME_AD_METH.EQ.2) CALL RK_DUCT_2            
          ELSEIF (NUM_PER_DIR.EQ.0) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CAV_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CAV_2            
          END IF
        END DO

C If you want to test the measurement extraction tool.
c	  call test_measurements
	  
        TIME=TIME+DELTA_T
	
	END DO

	write(6,*) 'Starting Backwards March.'
	DELTA_T = -abs(DELTA_T)


      DO TIME_STEP = TIME_STEP-1, LAST_TIME, -1

C	  Backward Time Evolution of State
        DO RK_STEP=1,3
          IF (NUM_PER_DIR.EQ.3) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_PER_1
            IF (TIME_AD_METH.EQ.2) CALL RK_PER_2            
          ELSEIF (NUM_PER_DIR.EQ.2) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CHAN_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CHAN_2            
          ELSEIF (NUM_PER_DIR.EQ.1) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_DUCT_1
            IF (TIME_AD_METH.EQ.2) CALL RK_DUCT_2            
          ELSEIF (NUM_PER_DIR.EQ.0) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CAV_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CAV_2            
          END IF
	  END DO
	  
c        DO RK_STEP=1,3
c          IF (NUM_PER_DIR.EQ.3) THEN
c		    write(6,*) 'Periodic Case Not Supported'
C	        CALL RK_ADJ_PER
c          ELSEIF (NUM_PER_DIR.EQ.2) THEN
c            CALL RK_ADJ_CHAN
c          ELSEIF (NUM_PER_DIR.EQ.1) THEN
c			write(6,*) 'Duct Case Not Supported'
C            CALL RK_ADJ_DUCT
c          ELSEIF (NUM_PER_DIR.EQ.0) THEN
c		    write(6,*) 'Cavity Case Not Supported'
C            CALL RK_ADJ_CAV
c          END IF
c        END DO
	  
        TIME=TIME+DELTA_T
        FIRST_TIME=.FALSE.
     	
      END DO
	write(6,*) 'Backwards March Complete'

! Calculate and display the runtime for the simulation
      CALL DATE_AND_TIME (VALUES=TIME_ARRAY)
      END_TIME=TIME_ARRAY(5)*3600+TIME_ARRAY(6)*60
     &   +TIME_ARRAY(7)+TIME_ARRAY(8)*0.001
      WRITE(*,*) 'Elapsed Time (sec): ',end_time-start_time
      WRITE(*,*) 'Seconds per Iteration: '
     &     ,(end_time-start_time)/N_TIME_STEPS

      TIME_STEP=TIME_STEP+1
      CALL SAVE_FLOW(.TRUE.)
      CALL SAVE_STATS(.TRUE.)
      
	return
	END



C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine REGULARIZATION
C This code is designed to find statistical averages of L2 norms of the
C error in a backwards march of the Navier Stokes Equations regularized
C with a Mized Time-Space Derivative.  The code DOES NOT FIND GRAMMER
C MISTAKES, like the fragment before this sentance.
	INCLUDE 'header_batch'
	INTEGER WIN_WIDTH, DELTA_WIN, MAXWIN
	INTEGER REG_COUNT, MAX_REG_COUNT
      INTEGER ind, INITIALIZE_TIME, LAST_TIME, i, k, j, m
	LOGICAL FLAG
      real*8     store_U1(0:200-1,0:NX,0:NZ,0:NY)
     &         , store_U2(0:200-1,0:NX,0:NZ,0:NY)
     &         , store_U3(0:200-1,0:NX,0:NZ,0:NY)
     &         , store_L2(0:200-1,0:NX,0:NZ,0:NY)
     	real*8 scaler1, scaler2(0:200-1), scaler3(0:200-1)
	INTEGER statistics, N_statistics, blow_up(0:200-1)
	complex*16 store_CU1(0:NX,0:NZ,0:NY+1)
     &         , store_CU2(0:NX,0:NZ,0:NY+1)
     &         , store_CU3(0:NX,0:NZ,0:NY+1)
	real*8 min_reg, max_reg, delta_reg, explode
	parameter (N_statistics = 50000, explode = 2.0)

	write(*,*)
	write(*,*)  '                       REGULARIZATION CODE'
	write(*,*)

	open(unit=85,file="basics.dat",status="unknown")
	open(unit=86,file="energy.dat",status="unknown")

C	Change INITIALIZE so that the appropriate estimation field is 
C	being read into name.f
      CALL INITIALIZE
	call init_march
      MAXWIN = 200
	DELTA_WIN = 1
	MAX_REG_COUNT = 0
	max_reg = 0.1
	min_reg = 0.0
	delta_reg = (max_reg-min_reg)/real(MAX_REG_COUNT+1)

      write(85,*) N_statistics
	write(85,*) MAXWIN
	write(85,*) DELTA_WIN
	write(85,*) MAX_REG_COUNT
	write(85,*) max_reg
	write(85,*) min_reg
	write(85,*) NKX
	write(85,*) TNKZ
	
	do i = 0, NKX+1
	  write(85,*) kx(i)
	end do
	
	do k = 0, TNKZ
	  write(85,*) kz(k)
	end do
	

! Initialize START_TIME for run timing
      CALL DATE_AND_TIME (VALUES=TIME_ARRAY)
      START_TIME=TIME_ARRAY(5)*3600+TIME_ARRAY(6)*60
     &   +TIME_ARRAY(7)+TIME_ARRAY(8)*0.001

C A flag to determine if we are considering the first time-step

c	write(85,*) '**********************************************'	
c	write(85,*) '************ Part 1: General Case ************'	
c	write(85,*) '**********************************************'	

	do REG_COUNT = 0, MAX_REG_COUNT
      reg_mixed = real(REG_COUNT)*delta_reg
      FIRST_TIME=.TRUE.

c	reg_mixed = 0.05

	blow_up = 0	
	scaler2(0:MAXWIN-1) = 0.0
	do m = 0,MAXWIN-1
	do i = 0,NKX
	  do k = 0,TNKZ
	    do j = 0,NY
	      store_L2(m,i,k,j) = 0.0
	    end do
	  end do
	end do
	end do
	
      INITIALIZE_TIME=3500
	GUSTS = .TRUE.
	DELTA_T = abs(DELTA_T)
	DO TIME_STEP = 1,INITIALIZE_TIME
C        WRITE(6,*) 'Now ending TIME_STEP = ',TIME_STEP
C	  Forward Time Evolution of State
        DO RK_STEP=1,3
          IF (NUM_PER_DIR.EQ.3) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_PER_1
            IF (TIME_AD_METH.EQ.2) CALL RK_PER_2            
          ELSEIF (NUM_PER_DIR.EQ.2) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CHAN_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CHAN_2            
          ELSEIF (NUM_PER_DIR.EQ.1) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_DUCT_1
            IF (TIME_AD_METH.EQ.2) CALL RK_DUCT_2            
          ELSEIF (NUM_PER_DIR.EQ.0) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CAV_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CAV_2            
          END IF
        END DO

        FIRST_TIME=.FALSE.

	END DO	

      INITIALIZE_TIME=50
	do statistics = 1,N_statistics
C INITIALIZE FLOW TO GET A SMOOTH TRAJECTORY AND POPULATE MEASUREMENTS
	GUSTS = .TRUE.
	DELTA_T = abs(DELTA_T)
      FLAG = .FALSE.

	DO TIME_STEP = 1,INITIALIZE_TIME
C        WRITE(6,*) 'Now ending TIME_STEP = ',TIME_STEP
C	  Forward Time Evolution of State
        DO RK_STEP=1,3
          IF (NUM_PER_DIR.EQ.3) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_PER_1
            IF (TIME_AD_METH.EQ.2) CALL RK_PER_2            
          ELSEIF (NUM_PER_DIR.EQ.2) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CHAN_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CHAN_2            
          ELSEIF (NUM_PER_DIR.EQ.1) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_DUCT_1
            IF (TIME_AD_METH.EQ.2) CALL RK_DUCT_2            
          ELSEIF (NUM_PER_DIR.EQ.0) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CAV_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CAV_2            
          END IF
        END DO

	END DO

	do i = 0,NKX
	do k = 0,TNKZ
	do j = 0,NY
	  store_CU1(i,k,j)= CU1(i,k,j)
	  store_CU2(i,k,j)= CU2(i,k,j)
	  store_CU3(i,k,j)= CU3(i,k,j)
      end do
	end do
	end do
	
c	write(6,*) 'Flow has been initialized!'
	GUSTS = .FALSE.
	LAST_TIME = TIME_STEP
	DO TIME_STEP = TIME_STEP, TIME_STEP+MAXWIN-1

C	  Compute the L2 norm in Physical Space:
	  call fft_xzy_to_physical(CU1,U1)
        call fft_xzy_to_physical(CU2,U2)
        call fft_xzy_to_physical(CU3,U3)

	  do i = 0,NX
	  do k = 0,NZ
	  do j = 0,NY
	    store_U1(TIME_STEP-LAST_TIME,i,k,j)= U1(i,k,j)
	    store_U2(TIME_STEP-LAST_TIME,i,k,j)= U2(i,k,j)
	    store_U3(TIME_STEP-LAST_TIME,i,k,j)= U3(i,k,j)
        end do
	  end do
	  end do
	  scaler1 = 0.0

	  do i = 0,NX
	    do k = 0,NZ
	      do j = 0,NY
	      scaler1 = scaler1 + (U1(i,k,j)*U1(i,k,j)
     &                        +  U2(i,k,j)*U2(i,k,j)
     &                        +  U3(i,k,j)*U3(i,k,j))*(LX/NX)*(LZ/NZ)
	      end do	    
	    end do
	  end do
      
C	  Save the L2 norm in vector to determine if explosion later.
	  scaler3(TIME_STEP-LAST_TIME) = sqrt(scaler1)
	
	  call fft_xzy_to_fourier(U1,CU1)
	  call fft_xzy_to_fourier(U2,CU2)
        call fft_xzy_to_fourier(U3,CU3)


C	  Forward Time Evolution of State
        DO RK_STEP=1,3
          IF (NUM_PER_DIR.EQ.3) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_PER_1
            IF (TIME_AD_METH.EQ.2) CALL RK_PER_2            
          ELSEIF (NUM_PER_DIR.EQ.2) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CHAN_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CHAN_2            
          ELSEIF (NUM_PER_DIR.EQ.1) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_DUCT_1
            IF (TIME_AD_METH.EQ.2) CALL RK_DUCT_2            
          ELSEIF (NUM_PER_DIR.EQ.0) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CAV_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CAV_2            
          END IF
        END DO

        TIME=TIME+DELTA_T
	END DO

c	write(6,*) 'Starting Backwards March.'
	DELTA_T = -abs(DELTA_T)

      DO TIME_STEP = TIME_STEP-1, LAST_TIME, -1

C	  Backward Time Evolution of State
        DO RK_STEP=1,3
          IF (NUM_PER_DIR.EQ.3) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_PER_1
            IF (TIME_AD_METH.EQ.2) CALL RK_PER_2            
          ELSEIF (NUM_PER_DIR.EQ.2) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CHAN_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CHAN_2            
          ELSEIF (NUM_PER_DIR.EQ.1) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_DUCT_1
            IF (TIME_AD_METH.EQ.2) CALL RK_DUCT_2            
          ELSEIF (NUM_PER_DIR.EQ.0) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CAV_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CAV_2            
          END IF
	  END DO
	  
	  call fft_xzy_to_physical(CU1,U1)
        call fft_xzy_to_physical(CU2,U2)
        call fft_xzy_to_physical(CU3,U3)
	
	  scaler1 = 0.0
	
	  do i = 0,NX
	    do k = 0,NZ
	      do j = 0,NY
	      scaler1 = scaler1 + (U1(i,k,j)*U1(i,k,j)
     &                        +  U2(i,k,j)*U2(i,k,j)
     &                        +  U3(i,k,j)*U3(i,k,j))*(LX/NX)*(LZ/NZ)
	      end do	    
	    end do
	  end do
	  
	  if (sqrt(scaler1).ge.explode) then
	    FLAG = .TRUE.
	    goto 6
	  end if

	  do i = 0,NX
	    do k = 0,NZ
	      do j = 0,NY
		  store_U1(TIME_STEP-LAST_TIME,i,k,j) = 
     &                     U1(i,k,j)-store_U1(TIME_STEP-LAST_TIME,i,k,j)
		  store_U2(TIME_STEP-LAST_TIME,i,k,j) = 
     &                     U2(i,k,j)-store_U2(TIME_STEP-LAST_TIME,i,k,j)
		  store_U3(TIME_STEP-LAST_TIME,i,k,j) = 
     &                     U3(i,k,j)-store_U3(TIME_STEP-LAST_TIME,i,k,j)
	      end do
	    end do
	  end do	
	
		
	  scaler1 = 0.0
	
	  do i = 0,NX
	    do k = 0,NZ
	      do j = 0,NY
	      scaler1 = scaler1  
     &    +  (store_U1(TIME_STEP-LAST_TIME,i,k,j)
     &              *store_U1(TIME_STEP-LAST_TIME,i,k,j)
     &    +   store_U2(TIME_STEP-LAST_TIME,i,k,j)
     &              *store_U2(TIME_STEP-LAST_TIME,i,k,j)
     &    +   store_U3(TIME_STEP-LAST_TIME,i,k,j)
     &             *store_U3(TIME_STEP-LAST_TIME,i,k,j))*(LX/NX)*(LZ/NZ)
	      end do	    
	    end do
	  end do

	  scaler2(TIME_STEP-LAST_TIME) = scaler2(TIME_STEP-LAST_TIME) 
     &                    + (sqrt(scaler1)/scaler3(TIME_STEP-LAST_TIME))
	
	  call fft_xzy_to_fourier(U1,CU1)
	  call fft_xzy_to_fourier(U2,CU2)
        call fft_xzy_to_fourier(U3,CU3)
	  
	do i=0,NKX
	  do k=0,TNKZ
	    do j = 0,NY
	    store_L2(TIME_STEP-LAST_TIME,i,k,j) = 
     & store_L2(TIME_STEP-LAST_TIME,i,k,j)+(CU1(i,k,j)-store_CU1(i,k,j))
     &                               *conjg(CU1(i,k,j)-store_CU1(i,k,j))
     &                                   + (CU2(i,k,j)-store_CU2(i,k,j))
     &                               *conjg(CU2(i,k,j)-store_CU2(i,k,j))
     &                                   + (CU3(i,k,j)-store_CU3(i,k,j))
     &                               *conjg(CU3(i,k,j)-store_CU3(i,k,j))
	    end do
	  end do
	end do

	TIME=TIME-DELTA_T

      END DO

c	write(6,*) 'Backwards March Complete'


    6 if (FLAG) then
	  do TIME_STEP = TIME_STEP-1, LAST_TIME, -1
	    blow_up(TIME_STEP-LAST_TIME) = blow_up(TIME_STEP-LAST_TIME)+1
	  end do
	end if
	
	do i = 0,NKX
	do k = 0,TNKZ
	do j = 0,NY
	  CU1(i,k,j) = store_CU1(i,k,j)
	  CU2(i,k,j) = store_CU2(i,k,j)
	  CU3(i,k,j) = store_CU3(i,k,j)
      end do
	end do
	end do
	
	call poisson_p_per
	
		
	write(*,*) 'Iteration ', statistics, ' of ', N_statistics
     & , ' iterations.'      
c     &    '( Regularization Constant = ', reg_mixed,
c     & ',    Window Width = ', WIN_WIDTH, ' )',
c     & '        Ave. Normalized L2 Err: ', scaler2/(statistics-blow_up),	
c     & '       Blow Up(s): ',blow_up	

C	END LOOP THROUGH STATISTICS
	end do
	
C	Write DATA TO FILE:
	do m = 0, MAXWIN-1
	  write(85,*) m-MAXWIN
	  write(85,*) scaler2(m)/(statistics-blow_up(m))
	  write(85,*) blow_up(m)
	  do i = 0,NKX
	    do k = 0,TNKZ
	      write(86,*) store_L2(m,i,k,0)/(statistics-blow_up(m))
	    end do
	  end do
	end do

C	END LOOP THROUGH REGULARIZATION CONSTANT
	end do

c	write(85,*) 
c	write(85,*) '***************************************'	
c	write(85,*) '************ Part 1: DONE! ************'	
c	write(85,*) '***************************************'	
	
! Calculate and display the runtime for the simulation
      CALL DATE_AND_TIME (VALUES=TIME_ARRAY)
      END_TIME=TIME_ARRAY(5)*3600+TIME_ARRAY(6)*60
     &   +TIME_ARRAY(7)+TIME_ARRAY(8)*0.001
      WRITE(*,*) 'Elapsed Time (sec): ',end_time-start_time
      WRITE(*,*) 'Seconds per Iteration: '
     &     ,(end_time-start_time)/(N_statistics*(MAX_REG_COUNT+1))

	close(85)
	close(86)
	
	return
	END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine cost_function(f_eval,alpha)
	real*8 f_eval,alpha,freq,amp, pi
	
	pi = 4.0*atan(1.0)
	
c	write(*,*) 'This is Cost Function!'
	amp = 100
	freq = 1000
	
	f_eval = amp*cos((pi/freq)*alpha)
	
	return
	end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine mnbrak(ax,bx,cx,fa,fb,fc)
C POINT BX is bracketed by AX and CX.
c Given a function func, and a given distinct initial points ax and bx, this
c routine searches in the downhill direction (defined by the function as 
c evaluated at the initial points) and returns new points ax, bx, and cx
c that bracket a minimum of the function.  Also returnd are the function
c values at the three points: fa, fb, and fc
c Parameters: GOLD is the default ration by which successive intervals
c   are magnified; GLIMIT is the maximum magnification allowed for a
c   parabolic-fit step.
c Algorithm adapted from Numerical Recipes in Fortran 77
c 2nd ED, pp. 393-394
c Adapted by Chris Colburn November 2007

	real*8 ax,bx,cx,fa,fb,fc,GOLD,GLIMIT,TINY
	parameter (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
	real*8 dum,fu,q,r,u,ulim
	
	call cost_function(fa,ax) ! fa=cost_function(ax)
	call cost_function(fb,bx) ! fb=cost_function(bx)
	if (fb.gt.fa) then
	  dum = ax
	  ax  = bx
	  bx  = dum
	  dum = fb
	  fb  = fa
	  fa  = dum
	end if
	
	cx = bx+GOLD*(bx-ax)
	call cost_function(fc,cx) ! fc = cost_function(cx)
	
	do while (fb.ge.fc)
	  r = (bx-ax)*(fb-fc)
	  q = (bx-cx)*(fb-fa)
	  u = bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
	  ulim = bx+GLIMIT*(cx-bx)
	  if ((bx-u)*(u-cx).gt.0.) then
	    call cost_function(fu,u) ! fu = cost_function(u)
	    if (fu.lt.fc) then
		ax = bx
		fa = fb
		bx = u
		fb = fu
		return
	    elseif (fu.gt.fb) then
		cx = u
		fc = fu
		return
	    end if
	    u = cx + GOLD*(cx-bx)
	    call cost_function(fu,u) ! fu = cost_function(u)
	  elseif ((cx-u)*(u-ulim).gt.0.) then
	    call cost_function(fu,u) ! fu = cost_function(u)
	    if (fu.lt.fc) then
		bx = cx
		cx = u
		u = cx+GOLD*(cx-bx)
		fb = fc
		fc = fu
		call cost_function(fu,u) ! fu = cost_function(u)
	    end if
	  elseif ((u-ulim)*(ulim-cx).gt.0.) then
	    u  = ulim
	    call cost_function(fu,u) ! fu = cost_function(u)
	  else
	    u = cx + GOLD*(cx-bx)
	    call cost_function(fu,u) ! fu = cost_function(u)
	  end if
	ax = bx
	bx = cx
	cx = u
	fa = fb
	fb = fc
	fc = fu
	end do
	
	return
	end


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine brent(ax,bx,cx,tol,xmin)
c Given a functon f, and given a bracketing triplet of abscissas ax,bx, cx
c (such that bx is between ax and cs and f(bx) is less than both f(ax) and
c f(cx)), this routine isolates the minimum to a fractional precision of
c about tol using Brent's method.  The abscissa of the minimum is returned
c as xmin, and the minimum function value is returned as breant.
c PARAMETERS: ITMAX is the maximum amount of iterations; GOLD is the golden
c   ratio; ZEPS is a small number that protects agains trying to achieve
c   fractional accuracy for a minimum that happens to be exactly zero.
c Algorithm adapted from Numerical Recipes in Fortran 77
c 2nd ED, pp. 397-398
c Adapted by Chris Colburn November 2007.

	integer iter, ITMAX
	real*8 ax,bx,cx,tol,xmin,CGOLD,ZEPS
	real*8 a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
	parameter (ITMAX=100, CGOLD=0.3819660, ZEPS=1.e-10)
	
	a=min(ax,cx)
	b=max(ax,cx)
	v=bx
	w=v
	x=v
	e=0.
	call cost_function(fx,x)    ! fx = f(x)
	fv=fx
	fw=fx
	
	do iter=1,ITMAX
	  xm = 0.5*(a+b)
	  tol1 = tol*abs(x)+ZEPS
	  tol2 = 2.*tol1
	  if(abs(x-xm).le.(tol2-0.5*(b-a))) goto 3
	  if(abs(e).gt.tol1) then
	    r=(x-w)*(fx-fv)
	    q=(x-v)*(fx-fw)
	    p=(x-v)*q-(x-w)*r
	    q=2.*(q-r)
	    if (q.gt.0.) p=-p
	    q=abs(q)
	    etemp=e
	    e=d
	    if(abs(p).ge.abs(0.5*q*etemp).or.p.le.q*(a-x).or.
     &						 p.ge.q*(b-x)) goto 1
	    d=p/q
	    u=x+d
	    if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
	    goto 2
	  end if
    1	  if (x.ge.xm) then
	    e=a-x
	  else
	    e=b-x
	  end if
	  d = CGOLD*e
    2	  if(abs(d).ge.tol1) then
	    u=x+d
	  else
	    u=x+sign(tol1,d)
	  end if
	  call cost_function(fu,u)	! fu=f(u)
	  if (fu.le.fx) then
	    if (u.ge.x) then
		a=x
	    else
		b=x
	    end if
	    v=w
	    fv=fw
	    w=x
	    fw=fx
	    x=u
	    fx=fu
	  else
	    if (u.lt.x) then
		a=u
	    else
		b=u
	    end if
	    if (fu.le.fw .or. w.eq.x) then
		v=w
		fv=fw
		w=u
		fw=fu
	    else if (fu.le.fv .or. v.eq.x .or. v.eq.w) then
		v=u
		fv=fu
	    end if
	  end if
	end do
	write(*,*) 'Brent Exceeded Maximum Iterations'
    3	xmin=x
    	write(*,*) iter
	
	return
	end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine test_penta
	include 'header_batch'
	
	integer i, NX1
	PARAMETER (NX1 = 0)
	real*8 dummy_a(0:NX1,0:ny),dummy_b(0:NX1,0:ny),dummy_c(0:NX1,0:ny)
	real*8 dummy_d(0:NX1,0:ny),dummy_e(0:NX1,0:ny),dummy_g(0:NX1,0:ny)
	real*8 store_a(0:NX1,0:ny),store_b(0:NX1,0:ny),store_c(0:NX1,0:ny)
	real*8 store_d(0:NX1,0:ny),store_e(0:NX1,0:ny)
	real*8 g_store(0:NX1,0:ny),test(0:NX1,0:ny)
	real*8 zbqlnor, mu, sigma, sigma2

	
	call zbqlini(0)
	mu = 0.0
	sigma = 15.0
	sigma2 = 2.*sigma
	
	if (ny.gt.4) then
				
	do i=0,NY
	  dummy_a(0,i) = ZBQLNOR(mu,sigma)
	  dummy_b(0,i) = ZBQLNOR(mu,sigma)
	  dummy_c(0,i) = ZBQLNOR(mu,sigma2)
	  dummy_d(0,i) = ZBQLNOR(mu,sigma)
	  dummy_e(0,i) = ZBQLNOR(mu,sigma)
	  dummy_g(0,i) = ZBQLNOR(mu,sigma)
	  if (i.le.1) dummy_a(0,i) = 0.
	  if (i.le.0) dummy_b(0,i) = 0.
	  if (i.ge.ny) dummy_d(0,i) = 0.
	  if (i.ge.ny-1) dummy_e(0,i) = 0.
	  store_a(0,i) = dummy_a(0,i)
	  store_b(0,i) = dummy_b(0,i)
	  store_c(0,i) = dummy_c(0,i)
	  store_d(0,i) = dummy_d(0,i)
	  store_e(0,i) = dummy_e(0,i)
	  g_store(0,i) = dummy_g(0,i)
C	write(*,'(I6,F15.4,F15.4,F15.4,F15.4,F15.4,F15.4)')i,dummy_a(0,i),
C     & dummy_b(0,i),dummy_c(0,i),dummy_d(0,i),dummy_e(0,i), dummy_g(0,i)
	end do
	
	call penta_real(dummy_a,dummy_b,dummy_c,dummy_d,dummy_e,dummy_g,
     & NY,NX1)
	
	i = 0
	test(0,0) = store_c(0,i)*dummy_g(0,i) +store_d(0,i)*dummy_g(0,i+1)
     &                                + store_e(0,i)*dummy_g(0,i+2) 
	write(*,'(I6,F15.4,F15.4,E15.4)') 0, test(0,0), g_store(0,0),
     & test(0,0)-g_store(0,0)
	
	i = 1
	test(0,1) = store_b(0,i)*dummy_g(0,i-1) +store_c(0,i)*dummy_g(0,i)
     &        + store_d(0,i)*dummy_g(0,i+1)+store_e(0,i)*dummy_g(0,i+2)
	write(*,'(I6,F15.4,F15.4,E15.4)') 1, test(0,1), g_store(0,1),
     & test(0,1)-g_store(0,1)
	
	do i=2,ny-2
	  test(0,i) = store_a(0,i)*dummy_g(0,i-2)+
     &              store_b(0,i)*dummy_g(0,i-1)+
     &              store_c(0,i)*dummy_g(0,i)  +
     &              store_d(0,i)*dummy_g(0,i+1)+
     &              store_e(0,i)*dummy_g(0,i+2)
	  write(*,'(I6,F15.4,F15.4,E15.4)') i, test(0,i), g_store(0,i),
     & test(0,i)-g_store(0,i)
	end do
	
	i = ny-1
	test(0,i) = store_a(0,i)*dummy_g(0,i-2)
     &              +store_b(0,i)*dummy_g(0,i-1)
     &              +store_c(0,i)*dummy_g(0,i)
     &              +store_d(0,i)*dummy_g(0,i+1)
     	write(*,'(I6,F15.4,F15.4,E15.4)') ny-1, test(0,ny-1), 
     & g_store(0,ny-1), test(0,ny-1)-g_store(0,ny-1)
	
	i = ny
	test(0,ny) = store_a(0,i)*dummy_g(0,i-2) 
     &           + store_b(0,i)*dummy_g(0,i-1)
     &           + store_c(0,i)*dummy_g(0,i)
     	write(*,'(I6,F15.4,F15.4,E15.4)') ny, test(0,ny), 
     & g_store(0,ny), test(0,ny)-g_store(0,ny)
     
      end if
     
      return
      end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine test_septa
	include 'header_batch'
	
	integer i, NX1
	PARAMETER (NX1 = 0)
	real*8 dummy_a(0:NX1,0:ny),dummy_b(0:NX1,0:ny),dummy_c(0:NX1,0:ny)
	real*8 dummy_d(0:NX1,0:ny),dummy_e(0:NX1,0:ny),dummy_f(0:NX1,0:ny)
	real*8 dummy_g(0:NX1,0:ny),dummy_h(0:NX1,0:ny)
	real*8 store_a(0:NX1,0:ny),store_b(0:NX1,0:ny),store_c(0:NX1,0:ny)
	real*8 store_d(0:NX1,0:ny),store_e(0:NX1,0:ny),store_f(0:NX1,0:ny)
	real*8 store_g(0:NX1,0:ny),store_h(0:NX1,0:ny),test(0:NX1,0:ny)
	real*8 zbqlnor, mu, sigma, sigma2

	
	call zbqlini(0)
	mu = 0.0
	sigma = 15.0
	sigma2 = 2.*sigma
	
	if (ny.gt.4) then
				
	do i=0,NY
	  dummy_a(0,i) = ZBQLNOR(mu,sigma)
	  dummy_b(0,i) = ZBQLNOR(mu,sigma)
	  dummy_c(0,i) = ZBQLNOR(mu,sigma)
	  dummy_d(0,i) = ZBQLNOR(mu,sigma2)
	  dummy_e(0,i) = ZBQLNOR(mu,sigma)
	  dummy_f(0,i) = ZBQLNOR(mu,sigma)
	  dummy_g(0,i) = ZBQLNOR(mu,sigma)
	  dummy_h(0,i) = ZBQLNOR(mu,sigma)
	  
	  if (i.le.2) dummy_a(0,i) = 0.
	  if (i.le.1) dummy_b(0,i) = 0.
	  if (i.le.0) dummy_c(0,i) = 0.
	  if (i.ge.ny) dummy_e(0,i) = 0.
	  if (i.ge.ny-1) dummy_f(0,i) = 0.
	  if (i.ge.ny-2) dummy_g(0,i) = 0.
	  
	  store_a(0,i) = dummy_a(0,i)
	  store_b(0,i) = dummy_b(0,i)
	  store_c(0,i) = dummy_c(0,i)
	  store_d(0,i) = dummy_d(0,i)
	  store_e(0,i) = dummy_e(0,i)
	  store_f(0,i) = dummy_f(0,i)
	  store_g(0,i) = dummy_g(0,i)
	  store_h(0,i) = dummy_h(0,i)
C	write(*,'(I6,F15.4,F15.4,F15.4,F15.4,F15.4,F15.4)')i,dummy_a(0,i),
C     & dummy_b(0,i),dummy_c(0,i),dummy_d(0,i),dummy_e(0,i), dummy_g(0,i)
	end do
	
	call septa_real(dummy_a,dummy_b,dummy_c,dummy_d,dummy_e,dummy_f,
     & dummy_g,dummy_h,NY,NX1)
	
	i = 0
	  test(0,i) = store_d(0,i)*dummy_h(0,i)   +
     &              store_e(0,i)*dummy_h(0,i+1) +
     &              store_f(0,i)*dummy_h(0,i+2) +
     &              store_g(0,i)*dummy_h(0,i+3)
	  write(*,'(I6,F15.4,F15.4,E15.4,F15.4)') i, test(0,i), 
     & store_h(0,i), test(0,i)-store_h(0,i), dummy_h(0,i)
	
	i = 1
	  test(0,i) = store_c(0,i)*dummy_h(0,i-1) +
     &              store_d(0,i)*dummy_h(0,i)   +
     &              store_e(0,i)*dummy_h(0,i+1) +
     &              store_f(0,i)*dummy_h(0,i+2) +
     &              store_g(0,i)*dummy_h(0,i+3)
	  write(*,'(I6,F15.4,F15.4,E15.4,F15.4)') i, test(0,i),
     & store_h(0,i), test(0,i)-store_h(0,i), dummy_h(0,i)
     
	i = 2
	  test(0,i) = store_b(0,i)*dummy_h(0,i-2) +
     &              store_c(0,i)*dummy_h(0,i-1) +
     &              store_d(0,i)*dummy_h(0,i)   +
     &              store_e(0,i)*dummy_h(0,i+1) +
     &              store_f(0,i)*dummy_h(0,i+2) +
     &              store_g(0,i)*dummy_h(0,i+3)
	  write(*,'(I6,F15.4,F15.4,E15.4,F15.4)') i, test(0,i),
     & store_h(0,i), test(0,i)-store_h(0,i), dummy_h(0,i)
     
	do i=3,ny-3
	  test(0,i) = store_a(0,i)*dummy_h(0,i-3) +
     &              store_b(0,i)*dummy_h(0,i-2) +
     &              store_c(0,i)*dummy_h(0,i-1) +
     &              store_d(0,i)*dummy_h(0,i)   +
     &              store_e(0,i)*dummy_h(0,i+1) +
     &              store_f(0,i)*dummy_h(0,i+2) +
     &              store_g(0,i)*dummy_h(0,i+3)
	  write(*,'(I6,F15.4,F15.4,E15.4,F15.4)') i, test(0,i),
     & store_h(0,i), test(0,i)-store_h(0,i), dummy_h(0,i)
	end do
	
	i = ny-2
	  test(0,i) = store_a(0,i)*dummy_h(0,i-3) +
     &              store_b(0,i)*dummy_h(0,i-2) +
     &              store_c(0,i)*dummy_h(0,i-1) +
     &              store_d(0,i)*dummy_h(0,i)   +
     &              store_e(0,i)*dummy_h(0,i+1) +
     &              store_f(0,i)*dummy_h(0,i+2)
	  write(*,'(I6,F15.4,F15.4,E15.4,F15.4)') i, test(0,i), 
     & store_h(0,i), test(0,i)-store_h(0,i), dummy_h(0,i)
	
	i = ny-1
	  test(0,i) = store_a(0,i)*dummy_h(0,i-3) +
     &              store_b(0,i)*dummy_h(0,i-2) +
     &              store_c(0,i)*dummy_h(0,i-1) +
     &              store_d(0,i)*dummy_h(0,i)   +
     &              store_e(0,i)*dummy_h(0,i+1)
	  write(*,'(I6,F15.4,F15.4,E15.4,F15.4)') i, test(0,i),
     & store_h(0,i), test(0,i)-store_h(0,i), dummy_h(0,i)
     

	i = ny
	  test(0,i) = store_a(0,i)*dummy_h(0,i-3) +
     &              store_b(0,i)*dummy_h(0,i-2) +
     &              store_c(0,i)*dummy_h(0,i-1) +
     &              store_d(0,i)*dummy_h(0,i)
	  write(*,'(I6,F15.4,F15.4,E15.4,F15.4)') i, test(0,i),
     & store_h(0,i), test(0,i)-store_h(0,i), dummy_h(0,i)
      end if
     
      return
      end




