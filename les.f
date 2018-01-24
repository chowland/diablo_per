      subroutine les_chan
C This subroutine models the terms owing to the subgrid scale stress
C if the computation is to be treated as an LES not a DNS
C This subroutine should be called when the velocity is in fourier space 
C in the periodic directions, on output, the velocity will be 
C in physical space.
C It is assumed that the test filter and the LES filter are performed
C by the same operation
C On output S1 should contain |S| which may be used again in les_chan_th
C if for the subgrid scalar dissipation

      include 'header_les'

      integer i,j,k,l,m,ij

      real*8 S1_mean(1:NY)
      real*8 NU_T_mean(1:NY)
 

! Variables for Dynamic Smagorinsky model:
      real*8 C_SMAG
      parameter (C_SMAG=0.13d0)
      real*8 DELTA_Y(0:NY+1),DELTA_YF(0:NY+1) 
      real*8 alpha,beta 
      real*8 denominator_sum

! Array to store the velocity index for each component of the strain rate tensor
      integer U_index1(6)
      integer U_index2(6)

! Here, alpha is the test/LES filter width ratio
      parameter (alpha=4.d0)
! beta is the LES/grid filter width ratio
      parameter (beta=3.d0)

! First, for all models, apply boundary conditions to the velocity field
! (fill ghost cells) to ensure accurate calculation of gradients
C Apply Boundary conditions to velocity field
      IF (USE_MPI) THEN
        CALL APPLY_BC_VEL_MPI
      ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
      END IF

      if (LES_MODEL_TYPE.EQ.1) then
!     Constant Smagorinsky model

! First, compute the rate of strain tensor S_ij

      call compute_strain_chan

! Compute |S| at GYF points, store in S1
! Interpolation to GYF points is easy since by definition
! GYF points are exactly midway between neighboring GY points
      DO J=1,NY
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=SQRT(
     &                2.d0*Sij(I,K,J,1)**2.d0
     &               +4.d0*(0.5d0*(Sij(I,K,J+1,4)+Sij(I,K,J,4)))**2.d0
     &               +4.d0*Sij(I,K,J,5)**2.d0
     &               +2.d0*Sij(I,K,J,2)**2.d0
     &               +4.d0*(0.5d0*(Sij(I,K,J+1,6)+Sij(I,K,J,6)))**2.d0 
     &               +2.d0*Sij(I,K,J,3)**2.d0 )
          END DO
        END DO
      END DO
! Extend |S| to ghost cells
      DO K=0,NZM
        DO I=0,NXM
          S1(I,K,0)=S1(I,K,1)
          S1(I,K,NY+1)=S1(I,K,NY)
        END DO
      END DO


! Now, compute |S|*S_ij, storing in Sij
! First compute at GYF points 
      DO J=1,NY
        DO K=0,NZM
          DO I=0,NXM
            Sij(I,K,J,1)=S1(I,K,J)*Sij(I,K,J,1)
            Sij(I,K,J,5)=S1(I,K,J)*Sij(I,K,J,5)
            Sij(I,K,J,2)=S1(I,K,J)*Sij(I,K,J,2)
            Sij(I,K,J,3)=S1(I,K,J)*Sij(I,K,J,3)
          END DO
        END DO
      END DO
! Now, compute at GY points, interpolating |S|
      DO J=1,NY+1
        DO K=0,NZM
          DO I=0,NXM
! |S| interpolated to GY point 
            TEMP(I,K,J)=(S1(I,K,J)*DYF(j-1)+S1(I,K,J-1)*DYF(j))
     &                  /(2.d0*DY(j))
! The terms dU1/dy and dU3/dy in CSij(:,:,:,4) and CSij(:,:,:,6) respectively
! are subtracted out from Sij here since they are treated implicitly
! in eddy viscosity terms
            Sij(I,K,J,4)=TEMP(I,K,J)*Sij(I,K,J,4)
            Sij(I,K,J,6)=TEMP(I,K,J)*Sij(I,K,J,6)
          END DO
        END DO
      END DO



! We now have |S|*S_ij stored in Sij in Physical space

! Convert |S|*S_ij to Fourier space
      do ij=1,6 
        CALL FFT_XZ_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0,NY+1)
      end do
! Compute the filter lengthscale
! Absorb -2.d0*C_SMAG**2.d0 here for effienciency
      DO J=1,NY
! At GYF points:
! Constant Smagorinsky
        DELTA_YF(J)=-2.d0*C_SMAG**2.d0
     &     *(DX(1)*beta*DYF(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
! Wall Damping
!        DELTA_YF(J)=
!     &    -2.d0*(0.1d0*(1.d0-exp((-GYF(J)/(NU*25.d0))**3.d0)))**2.d0
!     &            *(DX(1)*beta*DYF(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)

      END DO
! Extend to ghost cells 
      DELTA_YF(0)=DELTA_YF(1)
      DELTA_YF(NY+1)=DELTA_YF(NY)      

      DO J=1,NY+1
! At GY points:
! Constant Smagorinsky
        DELTA_Y(J)=-2.d0*C_SMAG**2.d0
     &        *(DX(1)*beta*DY(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
! Wall Damping
!        DELTA_Y(J)=
!     &    -2.d0*(0.1d0*(1.d0-exp((-GY(J)/(NU*25.d0))**3.d0)))**2.d0
!     &            *(DX(1)*beta*DY(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
      END DO

! Get the eddy viscosity at GY points
! NU_T = (C_S^2 * DELTA^2)*|S|

      DO J=1,NY+1
        DO K=0,NZM
          DO I=0,NXM
            NU_T(I,K,J)=-0.5d0*DELTA_Y(J)*TEMP(I,K,J)
          END DO
        END DO
      END DO

! Now, compute TAU, store in the corresponging Sij
      DO J=1,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CSij(I,K,J,1)=DELTA_YF(J)*CSij(I,K,J,1)
            CSij(I,K,J,5)=DELTA_YF(J)*CSij(I,K,J,5)
            CSij(I,K,J,2)=DELTA_YF(J)*CSij(I,K,J,2)
            CSij(I,K,J,3)=DELTA_YF(J)*CSij(I,K,J,3)
          END DO
        END DO
      END DO
      DO J=1,NY+1 
        DO K=0,TNKZ
          DO I=0,NKX
            CSij(I,K,J,4)=DELTA_Y(J)*CSij(I,K,J,4)
            CSij(I,K,J,6)=DELTA_Y(J)*CSij(I,K,J,6)
          END DO
        END DO
      END DO

! tau_ij is now contained in CSij in Fourier space

! Convert the velocity to physical space
      call FFT_XZ_TO_PHYSICAL(CU1,U1,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU2,U2,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU3,U3,0,NY+1)

      else if ((LES_MODEL_TYPE.EQ.2).or.(LES_MODEL_TYPE.eq.3)) then
! Here, use a dynamic smagorinsky model with or without scale similar part

! Compute the filter width
      DO J=1,NY
! At GYF points:
!        DELTA_YF(J)=(beta*DX(1)*DYF(J)*2.d0*beta*DZ(1))**(1.d0/3.d0)
        DELTA_YF(J)=sqrt((beta*DX(1))**2.d0+(DYF(J)*2.d0)**2.d0
     &      +(beta*DZ(1))**2.d0)
      END DO
! Extend to ghost cells
      DELTA_YF(0)=DELTA_YF(1)
      DELTA_YF(NY+1)=DELTA_YF(NY)

      DO J=1,NY+1
! At GY points:
!        DELTA_Y(J)=(beta*DX(1)*DY(J)*2.d0*beta*DZ(1))**(1.d0/3.d0)
       DELTA_Y(J)=sqrt((beta*DX(1))**2.d0+(DY(J)*2.d0)**2.d0
     &          +(beta*DZ(1))**2.d0)
      END DO

! We need to calculate the components of C, the dynamic coefficient
! C_DYN will be defined at GYF points

! Compute the rate of strain tensor, store in Sij
      call compute_strain_chan

 
! Compute |S| at GYF points, store in S1
! Interpolation to GYF points is easy since by definition
! GYF points are exactly midway between neighboring GY points
      DO J=1,min(J2+1,NY)
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=SQRT(
     &                2.d0*Sij(I,K,J,1)**2.d0
     &               +4.d0*(0.5d0*(Sij(I,K,J+1,4)+Sij(I,K,J,4)))**2.d0
     &               +4.d0*Sij(I,K,J,5)**2.d0
     &               +2.d0*Sij(I,K,J,2)**2.d0
     &               +4.d0*(0.5d0*(Sij(I,K,J+1,6)+Sij(I,K,J,6)))**2.d0 
     &               +2.d0*Sij(I,K,J,3)**2.d0 )
          END DO
        END DO
      END DO
! Extend |S| to ghost cells
      DO K=0,NZM
        DO I=0,NXM
          S1(I,K,0)=S1(I,K,1)
          S1(I,K,NY+1)=S1(I,K,NY)
        END DO
      END DO

! Convert Ui to physical space
      CALL FFT_XZ_TO_PHYSICAL(CU1,U1,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CU2,U2,0,NY+1)
      CALL FFT_XZ_TO_PHYSICAL(CU3,U3,0,NY+1)

      if (les_model_type.eq.3) then
! If we are using a dynmic mixed model, then calculate the
! filtered velocity to be used for the scale similar part
! Apply the filter to the LES velocity and save
      do j=0,NY+1
        do k=0,NZM
          do i=0,NXM
            U_2BAR(i,k,j,1)=U1(i,k,j)
            U_2BAR(i,k,j,3)=U3(i,k,j)
          end do
        end do
      end do
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
! Interpolate U2 to GYF points
            U_2BAR(i,k,j,2)=0.5d0*(U2(i,k,j+1)+U2(i,k,j))
          end do
        end do
      end do
      do k=0,NZM
        do i=0,NXM
! Extrapolate to GYF ghost cells
          U_2BAR(i,k,0,2)=U2(i,k,1)*(1.d0-DY(1)/(2.d0*DYF(1)))
     &                     -U2(i,k,2)*DY(1)/(2.d0*DYF(1))
          U_2BAR(i,k,NY+1,2)=U2(i,k,NY+1)*(1.d0+DY(NY+1)/(2.d0*DYF(NY)))
     &                     -U2(i,k,NY)*DY(NY+1)/(2.d0*DYF(NY)) 
        end do
      end do       
! Now, filter the velocity
      do i=1,3 
        call les_filter_chan(U_2BAR(0,0,0,i),0,NY+1)
      end do
      end if

! Compute C_DYN only every x # of timesteps
!      if (((MOD(TIME_STEP,10).eq.0).AND.(RK_STEP.eq.1))
!     &       .or.FIRST_TIME) THEN

! Filter |S| and store in S_2BAR
      DO J=0,NY+1 
        DO K=0,NZM
          DO I=0,NXM
            S_2BAR(I,K,J)=S1(I,K,J)
          END DO
        END DO
      END DO
      call les_filter_chan(S_2BAR,0,NY+1)

! Save a copy of the velocity which will be filtered twice more 

      if (les_model_type.eq.3) then 
! Do only if a scale similar part is needed
        do ij=1,3  
        do j=0,NY+1
          do k=0,NZ+1
            do i=0,NX+1
              U_4BAR(i,k,j,ij)=U_2BAR(i,k,j,ij)
            end do
          end do
        end do
        end do
        do i=1,3
          call les_filter_chan(U_4BAR(0,0,0,i),0,NY+1)
          call les_filter_chan(U_4BAR(0,0,0,i),0,NY+1)
        end do
      end if


! Zero C_DYN
      do j=0,NY+1
        C_DYN(j)=0.d0
      end do

! Set the velocity index for each component of the stress tensor
      U_index1(1)=1
      U_index2(1)=1

      U_index1(2)=2
      U_index2(2)=2 

      U_index1(3)=3
      U_index2(3)=3

      U_index1(4)=1
      U_index2(4)=2

      U_index1(5)=1
      U_index2(5)=3

      U_index1(6)=2
      U_index2(6)=3

! The prep. work is now done, we are ready to start the algorithm to compute
! the dynamic model coefficient

! Do over all non-repeating components of the stress tensor
      do ij=1,6

! Here        ij=1 -> l=1,m=1
!             ij=2 -> l=2,m=2
!             ij=3 -> l=3,m=3
!             ij=4 -> l=1,m=2
!             ij=5 -> l=1,m=3
!             ij=6 -> l=2,m=3       
  
! Zero the numerator and denominator:
        do j=0,NY+1
          do k=0,NZ+1
            do i=0,NX+1
              numerator(i,k,j)=0.d0
              denominator(i,k,j)=0.d0
            end do
          end do
        end do

! cross is used to multiply by two to include the contribution from the
! other symmetric term if we are dealing with a cross-term
        if (ij.le.3) then 
! We are computing a diagonal term
          cross=1.d0
        else
! We are computing a cross term
          cross=2.d0 
        end if 

! First, compute Mij

       if ((ij.eq.4).or.(ij.eq.6)) then
       do j=1,NY 
         do k=0,NZM
           do i=0,NXM
! Sij is defined at GY points, interpolate to GYF points           
               temp(i,k,j)=0.5d0*(Sij(i,k,j+1,ij)+Sij(i,k,j,ij))
           end do
         end do
       end do
       else
       do j=1,NY 
         do k=0,NZM
           do i=0,NXM
! Sij is defined at GYF points, no interpolation needed
             temp(i,k,j)=Sij(i,k,j,ij)
           end do
         end do
        end do
        end if
! Define temp at ghost points since it will be filtered
        do k=0,NZM
          do i=0,NXM
             temp(i,k,0)=temp(i,k,1)
             temp(i,k,NY+1)=temp(i,k,NY)
           end do
         end do
! Filter temp
       call les_filter_chan(temp,0,NY+1)
! Multiply by |S| filtered        
       do j=0,NY+1
         do k=0,NZM
           do i=0,NXM
             temp(i,k,j)=temp(i,k,j)*(alpha*DELTA_YF(j))**2.d0
     &                                    *S_2BAR(i,k,j) 
           end do
         end do
       end do
! Get second term of Mij
       if ((ij.eq.4).or.(ij.eq.6)) then
! Sij is defined at GY points, interpolate to GYF points
         do j=1,NY
           do k=0,NZM
             do i=0,NXM
               Mij(i,k,j)=DELTA_YF(j)**2.d0*S1(i,k,j)
     &               *0.5d0*(Sij(i,k,j+1,ij)+Sij(i,k,j,ij))
             end do
           end do
         end do
! Define Mij at ghost cells since it will be filtered
         do k=0,NZM
           do i=0,NXM
             Mij(i,k,0)=Mij(i,k,1)
             Mij(i,k,NY+1)=Mij(i,k,NY)
           end do
         end do
       else
! Sij is defined at GYF points, no interpolation needed
         do j=1,NY
           do i=0,NXM
             do k=0,NZM 
               Mij(i,k,j)=DELTA_YF(j)**2.d0*S1(i,k,j)*Sij(i,k,j,ij)
             end do
           end do
         end do
! Define Mij at ghost cells
         do i=0,NXM
           do k=0,NZM 
             Mij(i,k,0)=Mij(i,k,1)
             Mij(i,k,NY+1)=Mij(i,k,NY)
           end do
         end do
       end if

! Filter Mij
       call les_filter_chan(Mij,0,NY+1)
 
! Add the second term of Mij stored in temp
       do j=0,NY+1
         do k=0,NZM
           do i=0,NXM
             Mij(i,k,j)=temp(i,k,j)-Mij(i,k,j)    
           end do
         end do
       end do
! Now, compute Lij and add Lij*Mij to the numerator
! temp=Ui*Uj:
       SELECT CASE (ij)
       CASE(1)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U1(i,k,j)*U1(i,k,j)
             end do
           end do 
         end do
       CASE(2)
         do j=1,NY
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=0.25d0*(U2(i,k,j+1)+U2(i,k,j))**2.d0
             end do
           end do
         end do
! Get estimate at ghost cells, Use U2 at GY ghost cell directly
         do k=0,NZM
           do i=0,NXM
               temp(i,k,0)=U2(i,k,1)**2.d0
               temp(i,k,NY+1)=U2(i,k,NY+1)**2.d0
           end do 
         end do
       CASE(3)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U3(i,k,j)*U3(i,k,j)
             end do
           end do 
         end do
       CASE(4) 
         do j=1,NY
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U1(i,k,j)*0.5d0*(U2(i,k,j+1)+U2(i,k,j))
             end do
           end do
         end do
! Get estimate at ghost cells, Use U2 at GY ghost cell directly
         do k=0,NZM
           do i=0,NXM
             temp(i,k,0)=U1(i,k,0)*U2(i,k,1)
             temp(i,k,NY+1)=U1(i,k,NY+1)*U2(i,k,NY+1)
           end do 
         end do
       CASE(5)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U1(i,k,j)*U3(i,k,j)
             end do
           end do 
         end do
       CASE(6) 
         do j=1,NY
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=0.5d0*(U2(i,k,j+1)+U2(i,k,j))*U3(i,k,j)
             end do
           end do
         end do
! Get estimate at ghost cells, Use U2 at GY ghost cell directly
         do k=0,NZM
           do i=0,NXM
             temp(i,k,0)=U3(i,k,0)*U2(i,k,1)
             temp(i,k,NY+1)=U3(i,k,NY+1)*U2(i,k,NY+1)
           end do 
         end do
       END SELECT
! Filter temp
       call les_filter_chan(temp,0,NY+1)
! Add Lij*Mij to numerator
       do j=0,NY+1
         do k=0,NZM
           do i=0,NXM 
             numerator(i,k,j)=Mij(i,k,j)
     &        *(temp(i,k,j)
     &        -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij)))
           end do
         end do
       end do

       if (LES_MODEL_TYPE.eq.3) then
! If mixed model, include the scale similar part  
! Add Nij*Mij to the numerator piece-by-piece
! Term3 
         call les_filter_chan(temp,0,NY+1)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               numerator(i,k,j)=numerator(i,k,j)+temp(i,k,j)*Mij(i,k,j)
             end do
           end do
         end do
! Term 4
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U_2BAR(i,k,j,U_index1(ij))
     &                    *U_2BAR(i,k,j,U_index2(ij))
             end do
           end do
         end do
         call les_filter_chan(temp,0,NY+1)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
               numerator(i,k,j)=numerator(i,k,j)-temp(i,k,j)*Mij(i,k,j)
             end do
           end do
         end do 
! Terms 1 and 2
         call les_filter_chan(temp,0,NY+1)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM
              numerator(i,k,j)=numerator(i,k,j)
     &        -(temp(i,k,j)
     &         -U_4BAR(i,k,j,U_index1(ij))*U_4BAR(i,k,j,U_index2(ij)))
     &           *Mij(i,k,j)
             end do
           end do
         end do
       end if
! Now, the denominator for this ij
       do j=0,NY+1
         do k=0,NZM
           do i=0,NXM
             denominator(i,k,j)=Mij(i,k,j)*Mij(i,k,j)
           end do
         end do
       end do

! Get plane average of numerator and denominator, add to C
! Note, since both the numerator and denominator are averaged, there
! is not need to divide the sum by NX*NZ

       do j=0,min(NY+1,J2+1)
          denominator_sum=SUM(denominator(0:NXM,0:NZM,j))
          if (denominator_sum.ne.0.) then
            C_DYN(j)=C_DYN(j)
     &                 -0.5d0*cross*SUM(numerator(0:NXM,0:NZM,j))
     &                 /denominator_sum 
          else 
            C_DYN(j)=C_DYN(j)+0.d0  
          end if
       end do

! End to ij
      end do

! We are now done with the dynamic procedure to calculate C

! If C_DYN < 0 at any level, set C_DYN=0 for numerical stability
      do j=0,NY+1
        if (C_DYN(j).lt.0) C_DYN(j)=0.d0
      end do

! At this point we have C_DYN and Sij

! End if compute C_DYN
!       END IF

! Get the eddy viscosity at GY points
! NU_T = C_DYN * (DELTA^2)*|S|
! Use exact second order interpolation to get C_DYN and S1 interpolated to
! GY points
      DO J=J1,J2
        DO K=0,NZM
          DO I=0,NXM
            NU_T(I,K,J)=
     &        (DYF(j-1)*C_DYN(j)+DYF(j)*C_DYN(j-1))/(2.d0*DY(j))
     &        *DELTA_Y(j)**2.d0
     &        *(DYF(j-1)*S1(i,k,j)+DYF(j)*S1(i,k,j-1))/(2.d0*DY(j))
          END DO
        END DO
      END DO

! Calculate TAUij in physical space, stored in Sij

      do ij=1,6
      if (LES_MODEL_TYPE.eq.2) then
! Dynamic Smagorinsky model, no scale similar part
      if ((ij.eq.1).or.(ij.eq.2).or.(ij.eq.3).or.(ij.eq.5)) then
! Here, Sij is defined at GYF points 
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            Sij(i,k,j,ij)=
     &        -2.d0*C_DYN(j)*DELTA_YF(j)**2.d0*S1(i,k,j)*Sij(i,k,j,ij)
          end do 
        end do
      end do
      do k=0,NZM
        do i=0,NXM
! Get TAUij at ghost cells for computing derivative
          Sij(i,k,0,ij)=Sij(i,k,1,ij)
          Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
        end do
      end do
      else if ((ij.eq.4).or.(ij.eq.6)) then 
! Here, Sij is defined at GY points, interpolate C_DYN, etc 
! Use exact second order interpolation
! Sij(:,:,:,4)=0.5*(dU1/dy + dU2/dx)
! But dU1/dy will be accounted for as an implicit eddy viscosity term,
! So, subtract if off from Sij here
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            Sij(i,k,j,ij)=-2.d0
     &        *(DYF(j-1)*C_DYN(j)+DYF(j)*C_DYN(j-1))/(2.d0*DY(j))
     &        *DELTA_Y(j)**2.d0
     &        *(DYF(j-1)*S1(i,k,j)+DYF(j)*S1(i,k,j-1))/(2.d0*DY(j))
     &        *Sij(i,k,j,ij)
          end do
        end do
      end do 
! Get TAUij at ghost cells for computing derivative
      do k=0,NZM
        do i=0,NXM
          Sij(i,k,1,ij)=Sij(i,k,2,ij)
          Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
        end do
       end do

! End if ij
      end if

      else if (LES_MODEL_TYPE.eq.3) then
! Model type = 3, dynamic mixed model with scale similar part
! Always define temp at GYF points to match U_2BAR
! temp=Ui*Uj:
       SELECT CASE (ij)
       CASE(1)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM 
               temp(i,k,j)=U1(i,k,j)*U1(i,k,j)
             end do
           end do
         end do
       CASE(2)
         do j=1,NY
           do k=0,NZM
             do i=0,NXM 
! Interpolate U2^2 to GYF points
               temp(i,k,j)=(0.5d0*(U2(i,k,j)+U2(i,k,j+1)))**2.d0
             end do
           end do
         end do
! Get interpolate temp to ghost cells 
         do k=0,NZM
           do i=0,NXM 
             temp(i,k,0)=(U2(i,k,1)*(1.d0-DY(1)/(2.d0*DYF(1)))
     &                   -U2(i,k,2)*DY(1)/(2.d0*DYF(1)))**2.d0
             temp(i,k,NY+1)=(U2(i,k,NY+1)*(1.d0+DY(NY+1)
     &        /(2.d0*DYF(NY)))-U2(i,k,NY)*DY(NY+1)/(2.d0*DYF(NY)))**2.d0
           end do
         end do
       CASE(3)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM 
               temp(i,k,j)=U3(i,k,j)*U3(i,k,j)
             end do
           end do
         end do
       CASE(4) 
         do j=1,NY
           do k=0,NZM
             do i=0,NXM 
! Interpolate U2 to GYF points
               temp(i,k,j)=0.5d0*(U2(i,k,j)+U2(i,k,j+1))*U1(i,k,j)
             end do
           end do
         end do
! Get interpolate temp to ghost cells 
         do k=0,NZM
           do i=0,NXM 
             temp(i,k,0)=(U2(i,k,1)*(1.d0-DY(1)/(2.d0*DYF(1)))
     &                   -U2(i,k,2)*DY(1)/(2.d0*DYF(1)))*U1(i,k,0)
             temp(i,k,NY+1)=(U2(i,k,NY+1)*(1.d0+DY(NY+1)
     &        /(2.d0*DYF(NY)))-U2(i,k,NY)*DY(NY+1)/(2.d0*DYF(NY)))
     &         *U1(i,k,NY+1)
           end do
         end do
       CASE(5)
         do j=0,NY+1
           do k=0,NZM
             do i=0,NXM 
               temp(i,k,j)=U1(i,k,j)*U3(i,k,j)
             end do
           end do
         end do
       CASE(6) 
         do j=1,NY
           do k=0,NZM
             do i=0,NXM 
! Interpolate U3 to GY points
               temp(i,k,j)=0.5d0*(U2(i,k,j)+U2(i,k,j+1))*U3(i,k,j)
             end do
           end do
         end do
! Get interpolate temp to ghost cells 
         do k=0,NZM
           do i=0,NXM 
             temp(i,k,0)=(U2(i,k,1)*(1.d0-DY(1)/(2.d0*DYF(1)))
     &                  -U2(i,k,2)*DY(1)/(2.d0*DYF(1)))*U3(i,k,0)
             temp(i,k,NY+1)=(U2(i,k,NY+1)*(1.d0+DY(NY+1)
     &        /(2.d0*DYF(NY)))-U2(i,k,NY)*DY(NY+1)/(2.d0*DYF(NY)))
     &         *U3(i,k,NY+1)
           end do
         end do
       END SELECT


       call les_filter_chan(temp,0,NY+1)


      if ((ij.eq.1).or.(ij.eq.2).or.(ij.eq.3).or.(ij.eq.5)) then
! Here, Sij is defined at GYF points 
      do j=1,NY
        do k=0,NZM
          do i=0,NXM
            Sij(i,k,j,ij)=temp(i,k,j)
     &         -U_2BAR(i,k,j,U_index1(ij))*U_2BAR(i,k,j,U_index2(ij))
     &         -2.d0*C_DYN(j)*DELTA_YF(j)**2.d0*S1(i,k,j)*Sij(i,k,j,ij)
          end do
        end do 
      end do
! Get TAUij at ghost cells for computing derivative
      do k=0,NZM
        do i=0,NXM
          Sij(i,k,0,ij)=Sij(i,k,1,ij)
          Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
        end do
       end do
      else if ((ij.eq.4).or.(ij.eq.6)) then
! Here, Sij is defined at GY points, interpolate C_DYN, etc 
! Use exact second order interpolation
! Sij(:,:,:,4)=0.5*(dU1/dy + dU2/dx)
! But the dU1/dy term will be accounted for as an implicit eddy viscosity
! so subtract this term from Sij in the Smagorinsky part
      do j=2,NY
        do k=0,NZM
          do i=0,NXM
            Sij(i,k,j,ij)=
     &         (temp(i,k,j)*DYF(j-1)+temp(i,k,j-1)*DYF(j))/(2.d0*DY(j))
     &         -((U_2BAR(i,k,j,U_index1(ij))*DYF(j-1)
     &          +U_2BAR(i,k,j-1,U_index1(ij))*DYF(j))/(2.d0*DY(j)))
     &         *((U_2BAR(i,k,j,U_index2(ij))*DYF(j-1)
     &          +U_2BAR(i,k,j-1,U_index2(ij))*DYF(j))/(2.d0*DY(j)))
     &      -2.d0*(C_DYN(j)*DYF(j-1)+C_DYN(j-1)*DYF(j))/(2.d0*DY(j))
     &          *DELTA_Y(j)**2.d0
     &          *(S1(i,k,j)*DYF(j-1)+S1(i,k,j-1)*DYF(j))/(2.d0*DY(j))
     &          *Sij(i,k,j,ij)
          end do
        end do
      end do 
! Get TAUij at ghost cells for computing derivative
      do k=0,NZM
        do i=0,NXM
          Sij(i,k,1,ij)=Sij(i,k,2,ij)
          Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
        end do
       end do

! End if ij
       end if
! End if Mixed Model
       end if
! Convert TAUij, now stored in Sij to Fourier space
       call FFT_XZ_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0,NY+1)
! End do ij
       end do 

! End if LES_MODEL_TYPE dynamic Smagorinsky or Dynamic mixed model
      else
        pause 'Error, unsupported LES_MODEL_TYPE chosen'
      end if


! When using a Near-wall model, don't use LES at the wall
      IF ((U_BC_LOWER.EQ.3).or.(W_BC_LOWER.EQ.3)) then
         J1=2
      ELSE 
         J1=1
      END IF

      IF ((U_BC_UPPER.EQ.3).or.(W_BC_UPPER.EQ.3)) then
         J2=NY-1
      ELSE
         J2=NY
      END IF

! Zero Tau above J2+1 where it is not used
      DO ij=1,6
        DO J=J2+2,NY+1
          DO K=0,TNKZ
            DO I=0,NKX
              CSij(I,K,J,ij)=0.d0
            END DO
          END DO
        END DO
      END DO

! Now, add the subgrid scale forcing to CFi
! (This includes the subgrid scale stress as an explicit R-K term

      DO J=J1,J2
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CF1(I,K,J)
     &                -CIKX(I)*CSij(I,K,J,1)
     &              -(CSij(I,K,J+1,4)-CSij(I,K,J,4))/DYF(J)
     &                -CIKZ(K)*CSij(I,K,J,5)
            CF3(I,K,J)=CF3(I,K,J)
     &                -CIKX(I)*CSij(I,K,J,5)
     &              -(CSij(I,K,J+1,6)-CSij(I,K,J,6))/DYF(J)
     &                -CIKZ(K)*CSij(I,K,J,3)   
          END DO
        END DO
      END DO
      DO J=2,NY
        DO K=0,TNKZ
          DO I=0,NKX
           CF2(I,K,J)=CF2(I,K,J)
     &                -CIKX(I)*CSij(I,K,J,4)
     &                -(CSij(I,K,J,2)-CSij(I,K,J-1,2))/DY(j)
     &                -CIKZ(K)*CSij(I,K,J,6)
          END DO
        END DO
      END DO

! Add the first Crank-Nicolson term involving terms k-1
!      DO J=J1,J2
!        DO K=0,TNKZ
!          DO I=0,NKX 
!            CR1(I,K,J)=CR1(I,K,J)-TEMP1
!     &                  *(CSij(I,K,J+1,4)-CSij(I,K,J,4))/DYF(J)
!          END DO
!        END DO
!      END DO
!      DO J=J1,J2
!        DO K=0,TNKZ
!          DO I=0,NKX 
!            CR3(I,K,J)=CR3(I,K,J)-TEMP1
!     &                  *(CSij(I,K,J+1,6)-CSij(I,K,J,6))/DYF(J)
!          END DO
!        END DO
!      END DO

! Periodically, output mean quantities
      IF ((MOD(TIME_STEP,SAVE_STATS_INT).EQ.0).AND.(RK_STEP.EQ.1)) THEN
! Get plane averages
        do j=1,NY
          S1_mean(j)=0.d0
          NU_T_mean(j)=0.d0
          do i=0,NXM
          do k=0,NZM
            S1_mean(j)=S1_mean(j)+S1(I,K,J)
            NU_T_mean(j)=NU_T_mean(j)+NU_T(I,K,J)
          end do
          end do
          S1_mean(j)=S1_mean(j)/dble(NX*NZ)
          NU_T_mean(j)=NU_T_mean(j)/dble(NX*NZ)
        end do

        open(42,file='mean_les.txt',form='formatted',status='unknown')
        write(42,*) TIME_STEP,TIME,DELTA_T
        do j=1,NY
          write(42,420) j,GYF(J)
     &    ,(dble(CSij(0,0,J,ij)),ij=1,6),
!     &       C_DYN(J)*DELTA_YF(J)**2.d0*S1_mean(J),
     &       C_DYN(J),
     &      NU_T_mean(J)
        end do
      END IF
420     format(I3,' ',9(F20.9,' '))

      RETURN
      END



      subroutine compute_strain_chan
C This subroutine computes S_ij for the filtered velocity field
C The input velocity field should be in fourier space in the periodic
C directions.
C For use in the LES model in channel flow (2 periodic directions)
      include 'header_les'

      integer I,J,K,ij

       
      DO J=1,NY
        DO K=0,TNKZ
          DO I=0,NKX
            CSij(I,K,J,1)=CIKX(I)*CU1(I,K,J)
            CSij(I,K,J,2)=(CU2(I,K,J+1)-CU2(I,K,J))/DYF(j) 
            CSij(I,K,J,3)=CIKZ(K)*CU3(I,K,J)  
            CSij(I,K,J,5)=0.5d0*(CIKZ(K)*CU1(I,K,J)+CIKX(I)*CU3(I,K,J))
          END DO
        END DO
      END DO
      DO J=1,NY+1
        DO K=0,TNKZ
          DO I=0,NKX
            CSij(I,K,J,4)=0.5d0*( (CU1(I,K,J)-CU1(I,K,J-1))/DY(j)
     &                          +CIKX(I)*CU2(I,K,J) ) 
            CSij(I,K,J,6)=0.5d0*( CIKZ(K)*CU2(I,K,J)
     &                         +(CU3(I,K,J)-CU3(I,K,J-1))/DY(j) )
          END DO
        END DO
      END DO  


! Convert rate of strain tensor to physical space
      do ij=1,6
        call FFT_XZ_TO_PHYSICAL(CSij(0,0,0,ij),Sij(0,0,0,ij),0,NY+1)
      end do 
! We now have S_ij in Physical space
      RETURN
      END

 
      subroutine les_filter_chan(A,jstart,jend)
! This subroutine applies the les filter to the input field
! The indices to the start and end of the array in the y-direction
! are also inputted to make the routine cablable of filtering fields
! at either GYF or GY points.
! The array that is passed should be in physical space
      integer jstart,jend,NX,NY,NZ,NXM,NZM,N_TH
      INCLUDE 'grid_def'

      real*8 A(0:NX+1,0:NZ+1,0:NY+1)
      real*8 B(0:NX-1,0:NZ-1,0:NY+1)

      integer im2(0:NX-1),im1(0:NX-1),ip1(0:NX+1),ip2(0:NX+2)
      integer km2(0:NZ-1),km1(0:NZ-1),kp1(0:NZ+1),kp2(0:NZ+2)

! These are the weights for the filtering operation used
      real*8 W0,W1,W2,Wm1,Wm2,Wm1_j,W0_j,W1_j


! The following is for the 3-point trapezoidal rule, alpha*beta=sqrt(6)
!      Wm2=0.d0
!      Wm1=1.d0/4.d0
!      W0=1.d0/2.d0
!      W1=1.d0/4.d0
!      W2=0.d0
      Wm1_j=1.d0/4.d0  
      W0_j=1.d0/2.d0
      W1_j=1.d0/4.d0
! The following is for the 5-point trapezoidal rule, alpha*beta=9
      Wm2=1.d0/8.d0
      Wm1=1.d0/4.d0
      W0=1.d0/4.d0
      W1=1.d0/4.d0
      W2=1.d0/8.d0

      NXM=NX-1
      NZM=NZ-1

!      do j=0,NY+1
!        do k=0,NZM
!          do i=0,NXM
!            B(i,k,j)=A(i,k,j)
!          end do
!        end do
!      end do

! Filter in the periodic directions using cshift
! Note, cshift if not used since it appears to be slower
! Apply filter in the x-direction
!      B=Wm2*CSHIFT(B,-2,1)+Wm1*CSHIFT(B,-1,1)+W0*B+W1*CSHIFT(B,1,1)
!     &       +W2*CSHIFT(B,2,1)

! Filter using more efficient F77 syntax:
! Set up array to loop around periodic directions
      do i=2,NXM
        im2(i)=i-2
      end do
      im2(1)=NXM
      im2(0)=NX-2
      do i=1,NXM
        im1(i)=i-1
      end do
      im1(0)=NXM
      do i=0,NX-2
        ip1(i)=i+1
      end do
      ip1(NXM)=0
      do i=0,NX-3
        ip2(i)=i+2    
      end do
      ip2(NX-2)=0
      ip2(NXM)=1

      do j=jstart,jend
        do k=0,NZM
          do i=0,NXM
            B(i,k,j)=Wm2*A(im2(i),k,j)+Wm1*A(im1(i),k,j)+W0*A(i,k,j)
     &         +W1*A(ip1(i),k,j)+W2*A(ip2(i),k,j)
          end do
        end do  
      end do
 
! Apply filter in the z-diretion
!      B=Wm2*CSHIFT(B,-2,2)+Wm1*CSHIFT(B,-1,2)+W0*B+W1*CSHIFT(B,1,2)
!     &       +W2*CSHIFT(B,2,2)
! Filter using more efficient F77 syntax:
! Set up array to loop around periodic directions
      do k=2,NZM
        km2(k)=k-2
      end do
      km2(1)=NZM
      km2(0)=NZ-2
      do k=1,NZM
        km1(k)=k-1
      end do
      km1(0)=NZM
      do k=0,NZ-2
        kp1(k)=k+1
      end do
      kp1(NZM)=0
      do k=0,NZ-3
        kp2(k)=k+2    
      end do
      kp2(NZ-2)=0
      kp2(NZM)=1

      do j=jstart,jend
        do k=0,NZM
          do i=0,NXM
            A(i,k,j)=Wm2*B(i,km2(k),j)+Wm1*B(i,km1(k),j)+W0*B(i,k,j)
     &         +W1*B(i,kp1(k),j)+W2*B(i,kp2(k),j)
          end do
        end do  
      end do

! Apply filter in the vertical direction at all physical cells
! (filter is not applied to ghost cells, but the values of the ghost cells
! is used for averaging)
!      B(:,:,jstart+1:jend-1)=W0*B(:,:,jstart:jend-2)
!     &                      +W1*B(:,:,jstart+1:jend-1)
!     &                      +W2*B(:,:,jstart+2:jend)
! Use more efficient F77 syntax:
!       do j=jstart+1,jend-1
!         do k=0,NZM
!           do i=0,NXM
!             B(i,k,j)=Wm1_j*B(i,k,j-1)+W0_j*B(i,k,j)+W1_j*B(i,k,j+1)  
!           end do
!         end do
!       end do

!      do j=jstart,jend
!        do k=0,NZM
!          do i=0,NXM
!            A(i,k,j)=B(i,k,j)
!          end do
!        end do
!      end do

      return
      end

      subroutine les_chan_th(n)
C This subroutine models the subgridscale terms 
c in the scalar advection equation for scalar number n
C if the computation is to be treated as an LES not a DNS
C This subroutine should be called when the velocity is in fourier space 
C   in the periodic directions
C S1 should contain |S| which was calculated in les_chan

      include 'header_les'

      integer i,j,k,l,m,ij,n

      real*8 S1_mean(1:NY)
 

! Variables for Dynamic Smagorinsky model:
      real*8 C_SMAG
      parameter (C_SMAG=0.13d0)
      real*8 DELTA_Y(0:NY+1),DELTA_YF(0:NY+1) 
      real*8 alpha,beta 
      real*8 denominator_sum

! Here, alpha is the test/LES filter width ratio
      parameter (alpha=4.d0)
! beta is the LES/grid filter width ratio
      parameter (beta=3.d0)

! First, for all models, apply boundary conditions to the velocity field
! (fill ghost cells) to ensure accurate calculation of gradients
      IF (USE_MPI) THEN
        CALL APPLY_BC_VEL_MPI
      ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
      END IF

      if (LES_MODEL_TYPE.EQ.1) then
!     Constant Smagorinsky model

! First, compute the rate of strain tensor S_ij

      call compute_scalar_grad(n)

! Now, compute |S|*dTH/dx_i, storing in Sij
! First compute at GYF points 
      DO J=1,NY
        DO K=0,NZM
          DO I=0,NXM
            Sij(I,K,J,1)=S1(I,K,J)*Sij(I,K,J,1)
            Sij(I,K,J,3)=S1(I,K,J)*Sij(I,K,J,3)
          END DO
        END DO
      END DO
! Convert |S|*S_ij to Fourier space
      ij=1 
      CALL FFT_XZ_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0,NY+1)
      ij=3
      CALL FFT_XZ_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0,NY+1)

! Sij(:,:,:,2) is added through an implicit eddy viscosity
      DO J=1,NY+1
        DO K=0,TNKZ
          DO I=0,NKX
             CSij(I,K,J,2)=0.d0
          END DO
        END DO
      END DO

! We now have |S|*dTH/dx_i stored in Sij(:,:,:,1..3) in Physical space

! Compute the filter lengthscale
! Absorb -2.d0*C_SMAG**2.d0 here for effienciency
      DO J=1,NY
! At GYF points:
! Constant Smagorinsky
        DELTA_YF(J)=-C_SMAG**2.d0
     &     *(DX(1)*beta*DYF(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
! Wall Damping
!        DELTA_YF(J)=
!     &    -1.d0*(0.1d0*(1.d0-exp((-GYF(J)/(NU*25.d0))**3.d0)))**2.d0
!     &            *(DX(1)*beta*DYF(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)

      END DO
! Extend to ghost cells 
      DELTA_YF(0)=DELTA_YF(1)
      DELTA_YF(NY+1)=DELTA_YF(NY)      

      DO J=1,NY+1
! At GY points:
! Constant Smagorinsky
        DELTA_Y(J)=-C_SMAG**2.d0
     &        *(DX(1)*beta*DY(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
! Wall Damping
!        DELTA_Y(J)=
!     &    -1.d0*(0.1d0*(1.d0-exp((-GY(J)/(NU*25.d0))**3.d0)))**2.d0
!     &            *(DX(1)*beta*DY(J)*2.d0*DZ(1)*beta)**(2.d0/3.d0)
      END DO

! Get the eddy diffusivity at GY points
! NU_T = (C_S^2 * DELTA^2)*|S| With |S| interpolated to GY points
      DO J=J1,J2
        DO K=0,NZM
          DO I=0,NXM
            KAPPA_T(I,K,J,N)=-1.d0*DELTA_Y(J)
     &       *(S1(I,K,J)*DYF(J-1)+S1(I,K,J-1)*DYF(J))/(2.d0*DY(J))
          END DO
        END DO
      END DO

! Now, compute TAU, store in the corresponging Sij
      DO K=0,TNKZ
        DO I=0,NKX
          DO J=1,NY
            CSij(I,K,J,1)=DELTA_YF(J)*CSij(I,K,J,1)
            CSij(I,K,J,3)=DELTA_YF(J)*CSij(I,K,J,3)
          END DO
        END DO
      END DO

! tau_ij is now contained in CSij in Fourier space

! Convert the scalar to physical space
      call FFT_XZ_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N),0,NY+1)

      else if ((LES_MODEL_TYPE.EQ.2).or.(LES_MODEL_TYPE.eq.3)) then
! Here, use a dynamic smagorinsky model
! Note, there is no scale similar model for the scalar,
! so model type choice 2 and 3 are identical for the scalar equation

! Compute the filter width
      DO J=1,NY
! At GYF points:
!        DELTA_YF(J)=(beta*DX(1)*DYF(J)*2.d0*beta*DZ(1))**(1.d0/3.d0)
        DELTA_YF(J)=sqrt((beta*DX(1))**2.d0+(DYF(J)*2.d0)**2.d0
     &      +(beta*DZ(1))**2.d0)
      END DO
! Extend to ghost cells
      DELTA_YF(0)=DELTA_YF(1)
      DELTA_YF(NY+1)=DELTA_YF(NY)

      DO J=1,NY+1
! At GY points:
!        DELTA_Y(J)=beta*(beta*DX(1)*DY(J)*2.d0*beta*DZ(1))**(1.d0/3.d0)
       DELTA_Y(J)=sqrt((beta*DX(1))**2.d0+(DY(J)*2.d0)**2.d0
     &          +(beta*DZ(1))**2.d0)
      END DO

! We need to calculate the components of C, the dynamic coefficient
! C_DYN_TH will be defined at GYF points

! Compute the scalar gradient, store in Sij(:,:,:,1..3)
      call compute_scalar_grad(n)

! Convert the scalar to physical space
      CALL FFT_XZ_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N),0,NY+1)

! Compute C_DYN_TH only every x # of timesteps
      if ((MOD(TIME_STEP,10).eq.0).OR.FIRST_TIME) THEN

! Store TH in Sij(:,:,:,4) and apply the test filter
      do j=0,NY+1
        do k=0,NZM
          do i=0,NXM
            Sij(i,k,j,4)=TH(i,k,j,n)
          end do
        end do
      end do
      call les_filter_chan(Sij(0,0,0,4),0,NY+1)

! Zero C_DYN_TH
      do j=0,NY+1
        C_DYN_TH(j,N)=0.d0
      end do

! Do over all non-repeating components of the scalar gradient
      do ij=1,3

! Zero the numerator and denominator:
        do j=0,NY+1
          do k=0,NZ+1
            do i=0,NX+1
              numerator(i,k,j)=0.d0
              denominator(i,k,j)=0.d0
            end do
          end do
        end do
        
! First, compute Mij

       do k=0,NZM
         do i=0,NXM
           do j=1,min(NY,J2+1)
             if (ij.eq.2) then
! Sij is defined at GY points, interpolate to GYF points           
               temp(i,k,j)=0.5d0*(Sij(i,k,j+1,ij)+Sij(i,k,j,ij))
             else
! Sij is defined at GYF points, no interpolation needed
               temp(i,k,j)=Sij(i,k,j,ij)
             end if
! Define temp at ghost points since it will be filtered
             temp(i,k,0)=temp(i,k,1)
             temp(i,k,NY+1)=temp(i,k,NY)
           end do
         end do
       end do
! Filter temp
       call les_filter_chan(temp,0,min(NY+1,J2+1))
! Multiply by |S| filtered        
       do j=0,min(NY+1,J2+1)
         do k=0,NZM
           do i=0,NXM
             temp(i,k,j)=temp(i,k,j)*(alpha*DELTA_YF(j))**2.d0
     &                                    *S_2BAR(i,k,j) 
           end do
         end do
       end do
! Get second term of Mij
       if (ij.eq.2) then
! Sij is defined at GY points, interpolate to GYF points
         do k=0,NZM
           do i=0,NXM
             do j=1,min(NY,J2+1)
               Mij(i,k,j)=DELTA_YF(j)**2.d0*S1(i,k,j)
     &               *0.5d0*(Sij(i,k,j+1,ij)+Sij(i,k,j,ij))
             end do
! Define Mij at ghost cells since it will be filtered
             Mij(i,k,0)=Mij(i,k,1)
             Mij(i,k,NY+1)=Mij(i,k,NY)
           end do
         end do
       else
! Sij is defined at GYF points, no interpolation needed
         do i=0,NXM
           do k=0,NZM 
             do j=1,min(NY,J2+1)
               Mij(i,k,j)=DELTA_YF(j)**2.d0*S1(i,k,j)*Sij(i,k,j,ij)
             end do
! Define Mij at ghost cells
             Mij(i,k,0)=Mij(i,k,1)
             Mij(i,k,NY+1)=Mij(i,k,NY)
           end do
         end do
       end if

! Filter Mij
       call les_filter_chan(Mij,0,min(NY+1,J2+1))
 
! Add the second term of Mij stored in temp
       Mij=temp-Mij    
! Now, compute Lij and add Lij*Mij to the numerator
! temp=Ui*Uj:
       SELECT CASE (ij)
       CASE(1)
         do j=0,min(NY+1,J2+1)
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U1(i,k,j)*TH(i,k,j,n)
             end do
           end do 
         end do
       CASE(2)
         do j=1,min(NY,J2+1)
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=TH(i,k,j,n)*0.5d0*(U2(i,k,j+1)+U2(i,k,j))
             end do
           end do
         end do
! Get estimate at ghost cells, Use U2 at GY ghost cell directly
         do k=0,NZM
           do i=0,NXM
             temp(i,k,0)=TH(i,k,0,n)*U2(i,k,1)
             temp(i,k,NY+1)=TH(i,k,NY+1,n)*U2(i,k,NY+1)
           end do 
         end do
       CASE(3)
         do j=0,min(NY+1,J2+1)
           do k=0,NZM
             do i=0,NXM
               temp(i,k,j)=U3(i,k,j)*TH(i,k,j,n)
             end do
           end do 
         end do
       END SELECT
! Filter temp
       call les_filter_chan(temp,0,NY+1)
! Add Lij*Mij to numerator
! Recall that Sij(:,:,:,4) holds TH_2BAR
       do j=0,min(NY+1,J2+1)
         do k=0,NZM
           do i=0,NXM 
             numerator(i,k,j)=Mij(i,k,j)
     &              *(temp(i,k,j)-Sij(i,k,j,4)*U_2BAR(i,k,j,ij))
           end do
         end do
       end do

! Now, the denominator for this ij  
       do j=0,min(NY+1,J2+1)
         do k=0,NZM
           do i=0,NXM
             denominator(i,k,j)=Mij(i,k,j)*Mij(i,k,j)
           end do
         end do
       end do

! Get plane average of numerator and denominator, add to C
! Note, since both the numerator and denominator are averaged, there
! is not need to divide the sum by NX*NZ

       do j=0,min(NY+1,J2+1)
          denominator_sum=SUM(denominator(0:NXM,0:NZM,j))
          if (denominator_sum.ne.0.) then
            C_DYN_TH(j,n)=C_DYN_TH(j,n)
     &                 -0.5d0*SUM(numerator(0:NXM,0:NZM,j))
     &                 /denominator_sum 
          else 
            C_DYN_TH(j,n)=C_DYN_TH(j,n)+0.d0  
          end if
       end do

! End to ij
      end do

! We are now done with the dynamic procedure to calculate C

! If C_DYN_TH < 0 at any level, set C_DYN_TH=0 for numerical stability
      do j=0,NY+1
        if (C_DYN_TH(j,n).lt.0) C_DYN_TH(j,n)=0.d0
      end do

! End if compute C_DYN_TH
       END IF

! Get the eddy diffusivity at GY points
! KAPPA_T = C_DYN_TH * (DELTA^2)*|S|
! Use exact second order interpolation to get C_DYN_TH and S1 interpolated to
! GY points
      DO J=1,min(NY+1,J2+1)
        DO K=0,NZM
          DO I=0,NXM
            KAPPA_T(I,K,J,N)=
     &        (DYF(j-1)*C_DYN_TH(j,n)+DYF(j)*C_DYN_TH(j-1,n))
     &            /(2.d0*DY(j))
     &        *DELTA_Y(j)**2.d0
     &        *(DYF(j-1)*S1(i,k,j)+DYF(j)*S1(i,k,j-1))/(2.d0*DY(j))
          END DO
        END DO
      END DO


! At this point we have C_DYN_TH and dTH/dx_i (stored in Sij(:,:,:,1...3)
! Calculate lambda_i in physical space, stored in Sij(:,:,:,1..3)

      do ij=1,3,2
! Note, ij=2 is added to the
! equations implicitly through an eddy diffusivity
! Dynamic Smagorinsky model, no scale similar part
! Here, Sij is defined at GYF points for ij=1,3 and GY points for ij=2
      do k=0,NZM
        do i=0,NXM
          do j=1,min(NY,J2+1)
            Sij(i,k,j,ij)=
     &        -C_DYN_TH(j,n)*DELTA_YF(j)**2.d0*S1(i,k,j)*Sij(i,k,j,ij)
          end do 
! Get TAUij at ghost cells for computing derivative
          Sij(i,k,0,ij)=Sij(i,k,1,ij)
          Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
        end do
       end do
      end do
      ij=2
      do k=0,NZM
        do i=0,NXM
          do j=1,min(NY,J2+1)
            Sij(i,k,j,ij)=
     &     -(C_DYN_TH(j,n)*DYF(J-1)+C_DYN_TH(j-1,n)*DYF(J))/(2.d0*DY(j))
     &     *DELTA_Y(j)**2.d0
     &    *(S1(i,k,j)*DYF(J-1)+S1(i,k,j-1)*DYF(J))/(2.d0*DY(J))
     &       *Sij(i,k,j,ij)
          end do 
! Get TAUij at ghost cells for computing derivative
          Sij(i,k,0,ij)=Sij(i,k,1,ij)
          Sij(i,k,NY+1,ij)=Sij(i,k,NY,ij)
        end do
      end do 

      do ij=1,3
! Convert TAUij, now stored in Sij to Fourier space
        call FFT_XZ_TO_FOURIER(Sij(0,0,0,ij),CSij(0,0,0,ij),0
     &          ,min(NY+1,J2+1))
       end do  

! End if LES_MODEL_TYPE dynamic Smagorinsky or Dynamic model
      else
        pause 'Error, unsupported LES_MODEL_TYPE chosen'
      end if

! Now, add the subgrid scale forcing to CFi
! (This includes the subgrid scale stress as an explicit R-K term
! Include only CSij terms 1 and 3 since term 2 is accounted for
! as an implicit eddy diffusivity through KAPPA_T
      DO J=J1,J2
        DO K=0,TNKZ
          DO I=0,NKX
            CFTH(I,K,J,n)=CFTH(I,K,J,n)
     &                -CIKX(I)*CSij(I,K,J,1)
     &                -CIKZ(K)*CSij(I,K,J,3)
          END DO
       END DO
      END DO

      DO J=J1,J2
        DO K=0,TNKZ
          DO I=0,NKX
            CRTH(I,K,J,n)=CRTH(I,K,J,n)-TEMP1
     &         *(CSij(I,K,J+1,2)-CSij(I,K,J,2))/DYF(J)
          END DO
        END DO
      END DO

! Periodically, output mean quantities
      IF ((MOD(TIME_STEP,SAVE_STATS_INT).EQ.0).AND.(RK_STEP.EQ.1)) THEN
       open(43,file='mean_les_th.txt',form='formatted',status='unknown')
       write(43,*) TIME_STEP,TIME,DELTA_T
       do j=1,NY
          write(43,420) j,GYF(J)
     &    ,(dble(CSij(0,0,J,ij)),ij=1,3),
     &    SUM(KAPPA_T(0:NXM,0:NZM,J,N))/DBLE(NX*NZ),
     &       C_DYN_TH(J,n)*DELTA_YF(J)**2.d0
       end do
      END IF
420     format(I3,' ',5(F20.9,' '))

      RETURN
      END


      subroutine compute_scalar_grad(n)
C This subroutine computes dTH/dx_i for the filtered scalar field
C The input velocity field should be in fourier space in the periodic
C directions.
C For use in the LES model in channel flow (2 periodic directions)
C Store in Sij(:,:,:,1..3)
      include 'header_les'

      integer I,J,K,ij,n

      DO K=0,TNKZ
        DO I=0,NKX
          DO J=1,NY
            CSij(I,K,J,1)=CIKX(I)*CTH(I,K,J,n)
            CSij(I,K,J,3)=CIKZ(K)*CTH(I,K,J,n)
          END DO
          DO J=1,NY+1
! Define the vertical gradient at GY points
            CSij(I,K,J,2)=(CTH(I,K,J,n)-CTH(I,K,J-1,n))/DY(j)
          END DO
        END DO
      END DO

! Convert the scalar gradients to physical space
      do ij=1,3
        call FFT_XZ_TO_PHYSICAL(CSij(0,0,0,ij),Sij(0,0,0,ij),0,NY+1)
      end do

! We now have dTH/dx_i in Physical space
      RETURN
      END


      subroutine les_filter_chan_fourier(A,jstart,jend)
! This subroutine applies the les filter to the input field
! The filter is a spectral cutoff filter
! The indices to the start and end of the array in the y-direction
! are also inputted to make the routine cablable of filtering fields
! at either GYF or GY points.
! The array that is passed should be in physical space
      integer jstart,jend,NX,NZ,NY,NXM,NZM,i,j,k,N_TH

      INCLUDE 'grid_def'

      real*8 PI
      integer NKX,NKZ,TNKZ

      real*8 KX(0:NX/3),KZ(0:2*(NZ/3))  

      real*8 A(0:NX+1,0:NZ+1,0:NY+1)

      real alpha

      real*8 B(0:NX+1,0:NZ+1,0:NY+1)

      complex*16 CB(0:NX/2,0:NZ+1,0:NY+1)

      equivalence (B,CB)

! Set the ratio of filter scales
      parameter (alpha=2.d0)

      NXM=NX-1
      NZM=NZ-1


      PI = 4. * ATAN(1.0)

      LX=PI
      LZ=2.d0*PI
   
! Get the wavenumber vectors:
        NKX=NX/3
        DO I=0,NKX
          KX(I)=I*(2.*PI)/LX
        END DO

        NKZ=NZ/3
        TNKZ=NKZ*2
        DO K=0,NKZ
          KZ(K)=K*(2.*PI)/LZ
        END DO
        DO K=1,NKZ
          KZ(TNKZ+1-K)=-K*(2.*PI)/LZ
        END DO

       do j=0,NY+1
         do k=0,NZM
           do i=0,NXM
             B(i,k,j)=A(i,k,j)
           end do
         enddo
       end do 


! Convert to fourier space
      call fft_xz_to_fourier(B,CB,jstart,jend)

! Perform the filtering
      do j=jstart,jend
        do k=0,TNKZ
          do i=0,NKX
            if (sqrt(KX(I)**2.d0+KZ(K)**2.d0)
     &       .gt.sqrt(KX(NKX)**2.d0+KZ(NKZ)**2.d0)/alpha) then
                CB(i,k,j)=0.d0
            end if
          end do
        end do
      end do

! Now, convert back to physical space
      call fft_xz_to_physical(CB,B,jstart,jend)

       do j=jstart,jend 
         do k=0,NZM
           do i=0,NXM
             A(i,k,j)=B(i,k,j)
           end do
         end do
       end do 

      return
      end





