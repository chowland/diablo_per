! periodic.f, the fully-periodic-box solvers for diablo.           VERSION 2.5
!
! These solvers were written by Tom Bewley and John Taylor.
!******************************************************************************|

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine init_per
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

include 'header'

real*8 alpha
integer i,j,k

! Define variables for the Geophysical case
pi = 4.d0*atan(1.d0)
phi = 2.d0*pi*phi/360.d0
gamma = 2.d0*pi*gamma/360.d0
c_sin = cos(phi)*sin(gamma)/sin(phi)
c_cos = cos(phi)*cos(gamma)/sin(phi)

if (f_type == 3) then
    do j = 0,tnky
        do k = 0,tnkz_s
            do i = 0,nkx_s
                call random_number(alpha)
                f_phase(i,k,j) = alpha*2.d0*pi
            end do
        end do
    end do
end if

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine rk_per_1
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Main time-stepping algorithm for the fully periodic case
! INPUTS  (in Fourier space):  CUi, P, and (if k>1) CFi at (k-1)  (for i=1,2,3)
! OUTPUTS (in Fourier space):  CUi, P, and (if k<3) CFi at (k)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

include 'header'

real*8 temp1, temp2, temp3, temp4, temp5, temp6
integer i,j,k,n

! Start with data local in the X-direction

! Compute the RHS in Fourier Space, CRi.

! First, define the constants used for time-stepping
temp1 = nu*h_bar(rk_step)/2.d0
temp2 = beta_bar(rk_step)*h_bar(rk_step)
temp3 = zeta_bar(rk_step)*h_bar(rk_step)
temp4 = h_bar(rk_step)

do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
! Start with the explicit part of the Crank-Nicolson viscous term and
!  the pressure gradient treated with Explicit Euler:
            temp5 = 1 - temp1*(kx2_s(i)+ky2(j)+kz2_s(k))**beta

            cr1(i,k,j) = temp5*cu1(i,k,j) - temp4*(cikx_s(i)*cp(i,k,j))
            cr2(i,k,j) = temp5*cu2(i,k,j) - temp4*(ciky(j)*cp(i,k,j))
            cr3(i,k,j) = temp5*cu3(i,k,j) - temp4*(cikz_s(k)*cp(i,k,j))
        end do
    end do
end do

! For each scalar, start with the explict part of the Crank-Nicolson
! diffusive term for each scalar
do n = 1,n_th
    do j = 0,tnky_th
        do k = 0,tnkz_s_th
            do i = 0,nkx_s_th
                temp6 = 1-(temp1/pr(n))*(kx2_s_th(i) &
                            +ky2_th(j)+kz2_s_th(k))**beta &
                            -reaction(n)*temp1
                crth(i,k,j,n) = temp6*cth(i,k,j,n)
            end do
        end do
    end do
end do

if (rk_step  >  1) then
    do j = 0,tnky
        do k = 0,tnkz_s
            do i = 0,nkx_s
! Add the term: ZETA_BAR(RK_STEP)*R(U(RK_STEP-1))
                cr1(i,k,j) = cr1(i,k,j)+temp3*cf1(i,k,j)
                cr2(i,k,j) = cr2(i,k,j)+temp3*cf2(i,k,j)
                cr3(i,k,j) = cr3(i,k,j)+temp3*cf3(i,k,j)
            end do
        end do
    end do
    do n = 1,n_th
        do j = 0,tnky_th
            do k = 0,tnkz_s_th
                do i = 0,nkx_s_th
                    crth(i,k,j,n) = crth(i,k,j,n)+temp3*cfth(i,k,j,n)
                end do
            end do
        end do
    end do
end if


! If we are considering a linear background scalar gradient then add
! the term owing to advection of the background state.
! This allows us to consider a mean scalar gradient (ie stratification)
! even though the vertical boundary conditions are periodic.
! (In this case the passive scalar is a perturbation from a linear
! gradient. This gradient and the vertical domain size are used to
! make the passive scalar nondimensional, so here the nondimensional
! gradient is equal to one
do n = 1,n_th
    if (background_grad(n)) then
! If there is a background scalar gradient add advection term:
! Start by setting CFTH to zero:
        cfth(:,:,:,n) = 0.d0
! Now, add the velocity advection only to velocity modes
        do j = 0,tnky_th
            do k = 0,tnkz_s_th
                do i = 0,nkx_s_th
                    cfth(i,k,j,n) = -cu2(i,k,j)
                end do
            end do
        end do
    else
! Otherwise don't
        cfth(:,:,:,n) = 0.d0
    end if
end do

! Set the CFx variables to zero here so they can be used as working
! variables in USER_RHS_PER_FOURIER
cf1(:,:,:) = 0.d0
cf2(:,:,:) = 0.d0
cf3(:,:,:) = 0.d0

! Add the forcing terms CFx, CFTH using the subroutine in user_rhs.f
! Note that CUx, CTH are stored in Fourier space here
call user_rhs_per_fourier


! For Geophysical applications, add the Coriolis term
! Note, on an f-plane under the traditional approximation, C_SIN=C_COS=0
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cf1(i,k,j) = cf1(i,k,j)+i_ro_tau* &
                            (c_sin*cu2(i,k,j)+cu3(i,k,j))
            cf2(i,k,j) = cf2(i,k,j)+i_ro_tau* &
                            (-1.d0*c_cos*cu3(i,k,j)+c_sin*cu1(i,k,j))
            cf3(i,k,j) = cf3(i,k,j)+i_ro_tau* &
                            (c_cos*cu2(i,k,j)-cu1(i,k,j))
!           cf1(i,k,j) = i_ro_tau*cu3(i,k,j)
!           cf2(i,k,j) = 0.d0
!           cf3(i,k,j) = -1.d0*i_ro_tau*cu1(i,k,j)
        end do
    end do
end do


! Transform the scalar concentration to physical space
do n = 1,n_th
    call fft_xzy_mpi_to_physical_th(cth(0,0,0,n),th(0,0,0,n))
end do


! Compute the nonlinear terms for the passive scalar equation
! Do this before the nonlinear momentum terms to use Fi as a working
!  array before using it for the momentum equation.

do n = 1,n_th
    call fft_xzy_mpi_to_physical_interp(cu1,sth1)
    do j = 0,ny_s_th
        do k = 0,nz_s_th
            do i = 0,nxm_th
                sth1(i,k,j) = sth1(i,k,j)*th(i,k,j,n)
            end do
        end do
    end do
    call fft_xzy_mpi_to_fourier_th(sth1,csth1)
    do j = 0,tnky_th
        do k = 0,tnkz_s_th
            do i = 0,nkx_s_th
                cfth(i,k,j,n) = cfth(i,k,j,n)-cikx_s_th(i)*csth1(i,k,j)
            end do
        end do
    end do

    call fft_xzy_mpi_to_physical_interp(cu2,sth1)
    do j = 0,ny_s_th
        do k = 0,nz_s_th
            do i = 0,nxm_th
                sth1(i,k,j) = sth1(i,k,j)*th(i,k,j,n)
            end do
        end do
    end do
    call fft_xzy_mpi_to_fourier_th(sth1,csth1)
    do j = 0,tnky_th
        do k = 0,tnkz_s_th
            do i = 0,nkx_s_th
                cfth(i,k,j,n) = cfth(i,k,j,n)-ciky_th(j)*csth1(i,k,j)
            end do
        end do
    end do

    call fft_xzy_mpi_to_physical_interp(cu3,sth1)
    do j = 0,ny_s_th
        do k = 0,nz_s_th
            do i = 0,nxm_th
                sth1(i,k,j) = sth1(i,k,j)*th(i,k,j,n)
            end do
        end do
    end do
    call fft_xzy_mpi_to_fourier_th(sth1,csth1)
    do j = 0,tnky_th
        do k = 0,tnkz_s_th
            do i = 0,nkx_s_th
                cfth(i,k,j,n) = cfth(i,k,j,n)-cikz_s_th(k)*csth1(i,k,j)
            end do
        end do
    end do
end do

! Add R-K terms for the TH equation to the RHS
do n = 1,n_th
    do j = 0,tnky_th
        do k = 0,tnkz_s_th
            do i = 0,nkx_s_th
                crth(i,k,j,n) = crth(i,k,j,n)+temp2*cfth(i,k,j,n)
            end do
        end do
    end do
end do
! The RHS vector for the TH equation is now ready

! Inverse transform to physical space to compute the nonlinear terms
call fft_xzy_mpi_to_physical(cu1,u1)
call fft_xzy_mpi_to_physical(cu2,u2)
call fft_xzy_mpi_to_physical(cu3,u3)

! Compute the nonlinear terms for the momentum equations
do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            s1(i,k,j) = u1(i,k,j)*u1(i,k,j)
        end do
    end do
end do
call fft_xzy_mpi_to_fourier(s1,cs1)
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cf1(i,k,j) = cf1(i,k,j)-cikx_s(i)*cs1(i,k,j)
        end do
    end do
end do

do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            s1(i,k,j) = u1(i,k,j)*u2(i,k,j)
        end do
    end do
end do
call fft_xzy_mpi_to_fourier(s1,cs1)
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cf1(i,k,j) = cf1(i,k,j)-ciky(j)*cs1(i,k,j)
            cf2(i,k,j) = cf2(i,k,j)-cikx_s(i)*cs1(i,k,j)
        end do
    end do
end do

do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            s1(i,k,j) = u1(i,k,j)*u3(i,k,j)
        end do
    end do
end do
call fft_xzy_mpi_to_fourier(s1,cs1)
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cf1(i,k,j) = cf1(i,k,j)-cikz_s(k)*cs1(i,k,j)
            cf3(i,k,j) = cf3(i,k,j)-cikx_s(i)*cs1(i,k,j)
        end do
    end do
end do

do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            s1(i,k,j) = u2(i,k,j)*u2(i,k,j)
        end do
    end do
end do
call fft_xzy_mpi_to_fourier(s1,cs1)
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cf2(i,k,j) = cf2(i,k,j)-ciky(j)*cs1(i,k,j)
        end do
    end do
end do

do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            s1(i,k,j) = u2(i,k,j)*u3(i,k,j)
        end do
    end do
end do
call fft_xzy_mpi_to_fourier(s1,cs1)
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cf2(i,k,j) = cf2(i,k,j)-cikz_s(k)*cs1(i,k,j)
            cf3(i,k,j) = cf3(i,k,j)-ciky(j)*cs1(i,k,j)
        end do
    end do
end do

do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            s1(i,k,j) = u3(i,k,j)*u3(i,k,j)
        end do
    end do
end do
call fft_xzy_mpi_to_fourier(s1,cs1)
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cf3(i,k,j) = cf3(i,k,j)-cikz_s(k)*cs1(i,k,j)
        end do
    end do
end do
! Done with the computation of nonlinear terms

! If the scalar is active (RI_TAU NE 0), add the bouyancy forcing term
!  as explicit R-K
do n = 1,n_th
    if (ri_tau(n) /= 0.d0) then
! Fisrt, convert back to Fourier space
        call fft_xzy_mpi_to_fourier_th(th(0,0,0,n),cth(0,0,0,n))
! Add only to the velocity modes:
! ASSUME that the fine scale scalar structure does not affect the velocity
! Note, this needs to be modified before it can be used with varying grid sizes
        do j = 0,tnky
            do k = 0,tnkz_s
                do i = 0,nkx_s
                    cf2(i,k,j) = cf2(i,k,j)+ri_tau(n)*cth(i,k,j,n)
                end do
            end do
        end do
    end if
end do

! Now, add the R-K terms to the RHS
! Note, this has already been done for the scalar field TH
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cr1(i,k,j) = cr1(i,k,j)+temp2*cf1(i,k,j)
            cr2(i,k,j) = cr2(i,k,j)+temp2*cf2(i,k,j)
            cr3(i,k,j) = cr3(i,k,j)+temp2*cf3(i,k,j)
        end do
    end do
end do
! Computation of CRi complete.

! If Variable timestepping and done with one full R-K step, update
! DELTA_T based on the specified CFL number
! Note, this change will not take effect until the next timestep
! since the TEMP variables have already been defined
if ((variable_dt).and.(rk_step == 3) &
        .and.(mod(time_step,update_dt) == 0)) then
    call courant_mpi
end if

cu1 = (0.d0,0.d0)
cu2 = (0.d0,0.d0)
cu3 = (0.d0,0.d0)
cth = (0.d0,0.d0)

! Now solve the implicit system for the intermediate field.
! (In the fully-periodic case, this is easy!)
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            temp5 = 1.d0+temp1*(kx2_s(i)+ky2(j)+kz2_s(k))**beta
            cu1(i,k,j) = cr1(i,k,j)/temp5
            cu2(i,k,j) = cr2(i,k,j)/temp5
            cu3(i,k,j) = cr3(i,k,j)/temp5
        end do
    end do
end do
do n = 1,n_th
    do j = 0,tnky_th
        do k = 0,tnkz_s_th
            do i = 0,nkx_s_th
                temp6 = 1+(temp1/pr(n))*(kx2_s_th(i) &
                            +ky2_th(j)+kz2_s_th(k))**beta &
                            + reaction(n)*temp1
                cth(i,k,j,n) = crth(i,k,j,n)/temp6
            end do
        end do
    end do
end do
! First step of the Fractional Step algorithm complete.


! Begin second step of the Fractional Step algorithm, making u divergence free.

! Compute varphi, store in the variable CR1, and project velocity field
call rem_div_per

! Then, update the pressure with phi.
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cp(i,k,j) = cp(i,k,j)+cr1(i,k,j)/temp4
        end do
    end do
end do
! Second step of the Fractional Step algorithm complete.


return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine rem_div_per
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer i,j,k
real*8  temp5

! Compute phi, store in the variable CR1.
! Note the coefficient H_BAR is absorbed into phi.

cr1 = (0.d0,0.d0)

do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            temp5 = -(kx2_s(i)+ky2(j)+kz2_s(k)+eps)
            cr1(i,k,j) = (cikx_s(i)*cu1(i,k,j)+ciky(j)*cu2(i,k,j)+ &
                  cikz_s(k)*cu3(i,k,j))/temp5
        end do
    end do
! Then update CUi to make velocity field divergence-free.
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cu1(i,k,j) = cu1(i,k,j)-cikx_s(i)*cr1(i,k,j)
            cu2(i,k,j) = cu2(i,k,j)-ciky(j)*cr1(i,k,j)
            cu3(i,k,j) = cu3(i,k,j)-cikz_s(k)*cr1(i,k,j)
        end do
    end do
end do

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine poisson_p_per
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer i,j,k


do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cr1(i,k,j) = cikx_s(i)*cu1(i,k,j)
            cr2(i,k,j) = cikz_s(k)*cu3(i,k,j)
        end do
    end do
end do

call fft_xzy_mpi_to_physical(cr1,r1)
call fft_xzy_mpi_to_physical(cr2,r2)

do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            p(i,k,j) = r1(i,k,j)*r1(i,k,j)+r1(i,k,j)*r2(i,k,j)+ &
                            r2(i,k,j)*r2(i,k,j)
        end do
    end do
end do
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cf1(i,k,j) = ciky(j)*cu1(i,k,j)
            cf2(i,k,j) = cikx_s(i)*cu2(i,k,j)
            cf3(i,k,j) = cikz_s(k)*cu1(i,k,j)
            cr1(i,k,j) = cikx_s(i)*cu3(i,k,j)
            cr2(i,k,j) = cikz_s(k)*cu2(i,k,j)
            cr3(i,k,j) = ciky(j)*cu3(i,k,j)
        end do
    end do
end do
call fft_xzy_mpi_to_physical(cf1,f1)
call fft_xzy_mpi_to_physical(cf2,f2)
call fft_xzy_mpi_to_physical(cf3,f3)
call fft_xzy_mpi_to_physical(cr1,r1)
call fft_xzy_mpi_to_physical(cr2,r2)
call fft_xzy_mpi_to_physical(cr3,r3)
do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
        p(i,k,j) = 2*(p(i,k,j)+f1(i,k,j)*f2(i,k,j) &
                            +f3(i,k,j)*r1(i,k,j) &
                            +r2(i,k,j)*r3(i,k,j))
        end do
    end do
end do
call fft_xzy_mpi_to_fourier(p,cp)

do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cp(i,k,j) = cp(i,k,j)/(kx2_s(i)+ky2(j)+kz2_s(k)+eps)
        end do
    end do
end do

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine create_flow_per
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer i, j, k
real*8 rnum1,rnum2,rnum3
real*8 k0,k_mag, a0, theta0, phi0
character*60 fname
real*8 b3,f03,E3,sigma3,a3,alpha,h3,beta3,N03


! For an initial vortex, define the location of the centerline
real*8 xc(0:ny+1),zc(0:ny+1)

call mpi_barrier(mpi_comm_world, ierror)

if (flavor  ==  'Basic') then

    write(*,*) 'creating new flow from scratch.'

! Initialize random number generator
    call random_seed

    if (ic_type == 0) then
! Initizlize the flow using a Taylor-Green vortex
! Nondimensionalize with U0 and 1/kappa

        do j = 0,ny_s
            do k = 0,nz_s
                do i = 0,nxm
!            if (((j+int(rank/np_s)*(ny_s+1)) <= ny)
!     &          .and.((k+mod(rank,np_s)*(nz_s+1)) <= nz)) then
! Add a random phase
!            U1(I,K,J) = cos(2*pi*(GY_S(J))/LY)
!     &               *cos(2*pi*(GX(I))/LX)
!     &               *SIN(2*pi*(GZ_S(K))/LZ)
!            U2(I,K,J) = 0.d0
!            U3(I,K,J) = -cos(2*pi*(GY_S(J))/LY)
!     &               *sin(2*pi*(GX(I))/LX)
!     &               *COS(2*pi*(GZ_S(K))/LZ)
!            else
                    U1(i,k,j) = 0.d0
                    U2(i,k,j) = 0.d0
                    U3(i,k,j) = 0.d0
!            end if
                end do
            end do
        end do
    else if (ic_type == 1) then
! Start with an ideal vortex centered in the domain
        do j = 0,ny_s
!        xc(j) = lx/2.+(lx/5.)*sin(2*pi*gy(j)/ly)
!        zc(j) = lz/2.+(lz/5.)*sin(2*pi*gy(j)/ly)
            xc(j) = lx/2.
            zc(j) = lz/2.
            do k = 0,nz_s
                do i = 0,nxm
                    if ((gx(i)-xc(j))**2.+(gz_s(k)-zc(j))**2. > 0.1) then
! If we aren't too close to the vortex center
                        u1(i,k,j) = -1.d0*(gz_s(k)-zc(j)) &
                                    /((gx(i)-xc(j))**2.+(gz_s(k)-zc(j))**2.)
                        u3(i,k,j) = 1.d0*(gx(i)-xc(j)) &
                                    /((gx(i)-xc(j))**2.+(gz_s(k)-zc(j))**2.)
                        u2(i,k,j) = 0.d0
                    else
! Otherwise:
                        u1(i,k,j) = -1.d0*(gz_s(k)-zc(j))/0.1
                        u3(i,k,j) = 1.d0*(gx(i)-xc(j))/0.1
                        u2(i,k,j) = 0.d0
                    end if
                    call random_number(rnum1)
                    call random_number(rnum2)
                    call random_number(rnum3)
                    u1(i,k,j) = u1(i,k,j)+(rnum1-0.5)*kick
                    u2(i,k,j) = u2(i,k,j)+(rnum1-0.5)*kick
                    u3(i,k,j) = u3(i,k,j)+(rnum1-0.5)*kick
                end do
            end do
        end do
    else if (ic_type == 2) then
! Initilize with a plane wave
        cu1 = 0.d0
        cu2 = 0.d0
        cu3 = 0.d0
        if (rank == 0) then
            cu1(1,0,1) = 0.1/sqrt(2.0)/2.0
!            cu1(0,0,1) = 0.5
            cu2(1,0,1) = -0.1/sqrt(2.0)/2.0
            cu1(0,2,0)  =  0.5
        end if
    else if (ic_type == 3) then
! Initialize with a GM spectrum of internal waves
        b3 = 1300./1e-2/sqrt(Ri_tau(1)/nu)
        f03 = 7.3e-5*1e2*sqrt(Ri_tau(1))
        N03 = 5.2e-3*1e2*sqrt(Ri_tau(1))
        E3 = 6.3e-5
        sigma3 = 0.468
        a3 = sqrt(Ri_tau(1)*b3*E3*f03*sqrt(Ri_tau(1)-f03**2)) &
                /sqrt(2*pi*sigma3)
        h3 = 1/sqrt(2.)
        do j = 0,tnky
            do k = 0,tnkz_s
                do i = 0,nkx_s
                    if ( (kx2_s(i)+kz2_s(k)+ky2(j) <= 100.) .and. &
                            (ky(j) /= 0)) then
                        if (kx2_s(i)+kz2_s(k) == 0) then ! shear flow component
                            cs1(i,k,j) = sqrt(2*b3*E3*Ri_tau(1)/pi/sigma3/h3**2)
                                            *(ky2(j)+9*Ri_tau(1)*pi**2/b3**2/N03**2)**(-0.5)*(atan(h3*
                                                sqrt(Ri_tau(1)-f03**2)/f03/sqrt(h3**2+ky2(j))))**(0.5)
                            call random_number(beta3)
                            call random_number(alpha)
                            cu1(i,k,j) = sqrt(beta3)*cexp(cmplx(0,2.*pi*alpha))
                                            *cs1(i,k,j)
                            call random_number(alpha)
                            cu3(i,k,j) = sqrt(1-beta3)*cexp(cmplx(0,2.*pi*alpha))
                                            *cs1(i,k,j)
                        else
                            call random_number(alpha)
                            alpha = 2.*pi*alpha ! random phase of each wave
                            cs1(i,k,j) = a3*cexp(cmplx(0,alpha))*ky(j)/ &
                                            (sqrt(kx2_s(i)+kz2_s(k))* &
                                            sqrt(kx2_s(i)+kz2_s(k)+ky2(j)) &
                                            *(Ri_tau(1)*(kx2_s(i)+kz2_s(k))+f03**2*ky2(j)) &
                                            *(ky2(j)+pi**2*9*Ri_tau(1)/N03**2/b3**2))**(0.5)
                            cu1(i,k,j) = cs1(i,k,j)*ky(j)*kx_s(i)/sqrt( &
                                            (kx2_s(i)+kz2_s(k)+ky2(j))*(kx2_s(i)+kz2_s(k)))
                            cu2(i,k,j) = cs1(i,k,j)*sqrt(kx2_s(i)+kz2_s(k))/sqrt( &
                                            kx2_s(i)+kz2_s(k)+ky2(j))
                            cu3(i,k,j) = cs1(i,k,j)*ky(j)*kz_s(k)/sqrt( &
                                            (kx2_s(i)+kz2_s(k)+ky2(j))*(kx2_s(i)+kz2_s(k)))
                            cth(i,k,j,1) = cs1(i,k,j)*ci/sqrt(Ri_tau(1))
                        end if
                    end if
                end do
            end do
        end do
    else if (ic_type == 4) then
        ! Initialise with large shear + monochromatic IGW
        a0 = steep/2  ! steepness
        i = kx4       ! kx index
        j = ky4       ! ky index
        k = kz4       ! kz index
        cu1 = 0.d0
        cu2 = 0.d0
        cu3 = 0.d0
        if (rank == 0) then
            cu1(0,0,0) = u04    ! constant mean flow
            cu1(0,0,1) = ci      ! vertical shear at wavenumber 1
            cu1(i,k,j) = ci*a0*sqrt(Ri_tau(1))*kx_s(i) &
                            /sqrt(kx2_s(i)+kz2_s(k))/sqrt(kx2_s(i)+ky2(j)+kz2_s(k))
            cu2(i,k,j) = -ci*a0*sqrt(Ri_tau(1))/ky(j) &
                            *sqrt(kx2_s(i)+kz2_s(k))/sqrt(kx2_s(i)+ky2(j)+kz2_s(k))
            cu3(i,k,j) = ci*a0*sqrt(Ri_tau(1))*kz_s(k) &
                            /sqrt(kx2_s(i)+kz2_s(k))/sqrt(kx2_s(i)+ky2(j)+kz2_s(k))
        end if
    else if (ic_type == 5) then
        a0 = steep/2  ! steepness
        i = kx4       ! kx index
        j = ky4       ! ky index
        k = kz4       ! kz index
        cu1 = 0.d0
        cu2 = 0.d0
        cu3 = 0.d0
        if (rank == 0) then
            cu1(0,0,0) = u04    ! constant mean flow
            cu1(0,0,1) = -ci      ! vertical shear at wavenumber 1
            cu1(i,k,j) = ci*a0*sqrt(Ri_tau(1))*kx_s(i) &
                            /sqrt(kx2_s(i)+kz2_s(k))/sqrt(kx2_s(i)+ky2(j)+kz2_s(k))
            cu1(i,k,tnky+1-j) = ci*a0*sqrt(Ri_tau(1))*kx_s(i) &
                            /sqrt(kx2_s(i)+kz2_s(k)) &
                            /sqrt(kx2_s(i)+ky2(tnky+1-j)+kz2_s(k))
            cu2(i,k,j) = -ci*a0*sqrt(Ri_tau(1))/ky(j) &
                            *sqrt(kx2_s(i)+kz2_s(k))/sqrt(kx2_s(i)+ky2(j)+kz2_s(k))
            cu2(i,k,tnky+1-j) = -ci*a0*sqrt(Ri_tau(1))/ky(tnky+1-j) &
                            *sqrt(kx2_s(i)+kz2_s(k)) &
                            /sqrt(kx2_s(i)+ky2(tnky+1-j)+kz2_s(k))
            cu3(i,k,j) = ci*a0*sqrt(Ri_tau(1))*kz_s(k) &
                            /sqrt(kx2_s(i)+kz2_s(k))/sqrt(kx2_s(i)+ky2(j)+kz2_s(k))
            cu3(i,k,tnky+1-j) = ci*a0*sqrt(Ri_tau(1))*kz_s(k) &
                                /sqrt(kx2_s(i)+kz2_s(k)) &
                                /sqrt(kx2_s(i)+ky2(tnky+1-j)+kz2_s(k))
        end if
! Initialise with large shear + standing IGW
!       do j = 0,NY_S
!         do k = 0,NZ_S
!           do i = 0,NXM
!             U1(i,k,j) = sin(GY_S(j))+sqrt(RI_TAU(1))
!  &                  /2/sqrt(5.0)*cos(4*GX(i))*sin(8*GY_S(j))
!             U2(i,k,j) = -sqrt(RI_TAU(1))/4/sqrt(5.0)
!  &                          *sin(4*GX(i))*cos(8*GY_S(j))
!           end do
!         end do
!       end do
    else
        write(*,*) 'warning, undefined ics in periodic.f'
    end if

    ! if ((ic_type < 2) .or. (ic_type == 5)) then
    if (ic_type < 2) then
        call fft_xzy_mpi_to_fourier(u1,cu1)
        call fft_xzy_mpi_to_fourier(u3,cu3)
        call fft_xzy_mpi_to_fourier(u2,cu2)
    end if

    if (rank == 0) then
        do j = 0,tnky
            do k = 0,tnkz_s
                do i = 0,nkx_s
                    call random_number(rnum1)
                    call random_number(rnum2)
                    rnum1 = (rnum1-0.5d0)*2.d0
                    rnum2 = (rnum2-0.5d0)*2.d0
                    cu1(i,k,j) = cu1(i,k,j)+cmplx(rnum1,rnum2)*kick
                    call random_number(rnum1)
                    call random_number(rnum2)
                    rnum1 = (rnum1-0.5d0)*2.d0
                    rnum2 = (rnum2-0.5d0)*2.d0
                    cu2(i,k,j) = cu2(i,k,j)+cmplx(rnum1,rnum2)*kick
                    call random_number(rnum1)
                    call random_number(rnum2)
                    rnum1 = (rnum1-0.5d0)*2.d0
                    rnum2 = (rnum2-0.5d0)*2.d0
                    cu3(i,k,j) = cu3(i,k,j)+cmplx(rnum1,rnum2)*kick
                    call random_number(rnum1)
                    call random_number(rnum2)
                    rnum1 = (rnum1-0.5d0)*2.d0
                    rnum2 = (rnum2-0.5d0)*2.d0
                end do
            end do
        end do
    end if

else
    write(*,*) 'Unknown flavour, flow-field not created'
end if

call mpi_barrier(mpi_comm_world,ierror)

if (rank == 0) write(*,*) 'Calling rem_div...'
call rem_div_per

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine create_th_per
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer i,j,k,n
real*8 a0,theta0,phi0

! Note, Since stratification is not permitted in the periodic flow field
! Any background stratification must be added to the governing equations

if ((ic_type == 0).or.(ic_type == 1).or.(ic_type == 3)) then
    do n = 1,n_th
        if (create_new_th(n)) then
            do j = 0,ny_s_th
                do k = 0,nz_s_th
                    do i = 0,nxm_th
                        th(i,k,j,n) = 0.d0
                    end do
                end do
            end do
        end if
    end do
else if (ic_type == 2) then
! for plane wave
    do n = 1,n_th
        if (create_new_th(n)) then
            cth(:,:,:,n) = 0.d0
            if ((rank == 0) .and. (n == 1)) then
                cth(1,0,1,n) = ci*0.1/2.0
            end if
        end if
    end do
else if (ic_type == 4) then
    a0 = steep/2  ! steepness
    i = kx4       ! kx index
    j = ky4       ! ky index
    k = kz4       ! kz index
    if (create_new_th(1) .and. rank == 0) then
        cth(i,k,j,1) = -a0/ky(j)
    end if
else if (ic_type == 5) then
    a0 = steep/2  ! steepness
    i = kx4       ! kx index
    j = ky4       ! ky index
    k = kz4       ! kz index
    if (create_new_th(1) .and. rank == 0) then
        cth(i,k,j,1) = -a0/ky(j)
        cth(i,k,tnky+1-j,1) = -a0/ky(tnky+1-j)
    end if
! if (create_new_th(1)) then
    ! do j = 0,ny_s_th
    !   do k = 0,nz_s_th
    !     do i = 0,nxm_th
    !       th(i,k,j,1) = 1/4.0*cos(4*gx(i))*cos(8*gy_s(j))
    !     end do
    !   end do
    ! end do
! end if
else
    write(*,*) 'Unknown ic_type in create_th_per'
end if

! Transfer to Fourier space
! if (ic_type /= 2 .and. ic_type /= 4) then
if ((ic_type /= 2) .and. (ic_type < 4)) then
    do n = 1,n_th
        call fft_xzy_mpi_to_fourier_th(th(0,0,0,n),cth(0,0,0,n))
    end do
end if

return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine create_grid_per
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer i,j,k

do i = 0,nx
    gx(i) = (i*lx)/nx
    dx(i) = lx/nx
    if (verbosity > 3) write(*,*) 'gx(',i,') = ',gx(i)
end do

do k = 0,nz
    gz(k) = (k*lz)/nz
    dz(k) = lz/nz
    if (verbosity > 3) write(*,*) 'gz(',k,') = ',gz(k)
end do

do j = 0,ny
    gy(j) = (j*ly)/ny
    dy(j) = ly/ny
    if (verbosity > 3) write(*,*) 'gy(',j,') = ',gy(j)
end do

do j = 0,ny_s
    if ((j+int(rank/np_s)*(ny_s+1)) <= ny) then
        gy_s(j) = gy(j+int(rank/np_s)*(ny_s+1))
        dy_s(j) = dy(j+int(rank/np_s)*(ny_s+1))
        if (verbosity > 3) write(*,*) 'gy_s(',j,') = ',gy_s(j)
    end if
end do
do k = 0,nz_s
    if ((k+mod(rank,np_s)*(nz_s+1)) <= nz) then
        gz_s(k) = gz(k+mod(rank,np_s)*(nz_s+1))
        dz_s(k) = dz(k+mod(rank,np_s)*(nz_s+1))
        if (verbosity > 3) write(*,*) 'gz_s(',k,') = ',gz_s(k)
    end if
end do


do i = 0,nx_th
    gx_th(i) = (i*lx)/nx_th
    dx_th(i) = lx/nx_th
    if (verbosity > 3) write(*,*) 'gx_th(',i,') = ',gx_th(i)
end do

do k = 0,nz_th
    gz_th(k) = (k*lz)/nz_th
    dz_th(k) = lz/nz_th
    if (verbosity > 3) write(*,*) 'gz_th(',k,') = ',gz_th(k)
end do

do j = 0,ny_th
    gy_th(j) = (j*ly)/ny_th
    dy_th(j) = ly/ny_th
    if (verbosity > 3) write(*,*) 'gy_th(',j,') = ',gy_th(j)
end do

do j = 0,ny_s_th
    if ((j+int(rank/np_s)*(ny_s_th+1)) <= ny_th) then
        gy_s_th(j) = gy_th(j+int(rank/np_s)*(ny_s_th+1))
        dy_s_th(j) = dy_th(j+int(rank/np_s)*(ny_s_th+1))
    end if
end do
do k = 0,nz_s_th
    if ((k+mod(rank,np_s)*(nz_s_th+1)) <= nz_th) then
        gz_s_th(k) = gz_th(k+mod(rank,np_s)*(nz_s_th+1))
        dz_s_th(k) = dz_th(k+mod(rank,np_s)*(nz_s_th+1))
    end if
end do

return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine input_per
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
real    version, current_version
integer n

! Read in input parameters specific for periodic flow case
open (11, file = 'input_per.dat', form = 'formatted', status = 'old')
! Read input file.

current_version = 2.5
read(11,*)
read(11,*)
read(11,*)
read(11,*)
read(11,*) version
if (version  /=  current_version) &
    stop 'wrong input data format input_per'
if (rank == 0) write(*,*) 'version: ',version
read(11,*)
read(11,*) ic_type, kick
if (rank == 0) write(*,*) 'ic_type,kick: ',ic_type,kick
read(11,*)
read(11,*) i_Ro_tau, phi, gamma, g_tau, beta
if (rank == 0) write(*,*) 'i_Ro_tau,phi,gamma,g_tau,beta: ' &
    ,i_Ro_tau,phi,gamma,g_tau,beta
read(11,*)
read(11,*) f_type, force_shear, target_reb
do n = 1,n_th
    read(11,*)
    read(11,*) background_grad(n)
    if (rank == 0) write(*,*) 'background_grad(n): ', background_grad(n)
    read(11,*)
    read(11,*) u04, steep, kx4, ky4, kz4
end do

return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine save_stats_per(final)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
character*60 fname
logical final

integer i,j,k,n, iz, iy
real*8 uc, ubulk

integer irank, ranky_mov, rankz_mov

real*8 ubar,vbar,wbar
real*8 thbar(1:n_th)

real*8 u1rms,u2rms,u3rms
real*8 u1rms_sum,u2rms_sum,u3rms_sum
real*8 thrms_sum(1:n_th)
real*8 l_th(1:n_th),l_th_sum(1:n_th),thflux_sum(1:n_th)

real*8 u1me(0:ny_s), u3me(0:ny_s), th1me(0:ny_s)
real*8 u1u2(0:ny_s), u3u2(0:ny_s), thu2(0:ny_s_th,1:n_th)
real*8 u1u1(0:ny_s), u2u2(0:ny_s), u3u3(0:ny_s)
real*8 u1rms_h(0:ny_s), u2rms_h(0:ny_s), u3rms_h(0:ny_s)
real*8 thth(0:ny_s_th,1:n_th), thth_h(0:ny_s_th,1:n_th)
real*8 u1u2_sum(0:ny_s), u3u2_sum(0:ny_s)
real*8 thu2_sum(0:ny_s_th,1:n_th)
real*8 e_l(0:tnky), e_s(0:tnky), spectrum(0:tnky), e(0:tnky)
character(10) :: gname

call mpi_barrier(mpi_comm_world,ierror)

if (rank == 0) then
    write(*,*) 'Time , delta_t = ', time, delta_t
end if

if ((flavor == 'Basic').or.(flavor == 'Chemotaxis')) then

! Note that this routine uses CRi and CRTH for storage, so it should
! only be used between full R-K timesteps

    if (rank == 0) then
        write(*,*) 'Saving flow statistics.'
    end if

! Store the velocity in Fourier space in CRi
    do j = 0,tnky
        do k = 0,tnkz_s
            do i = 0,nkx_s
                cr1(i,k,j) = cu1(i,k,j)
                cr2(i,k,j) = cu2(i,k,j)
                cr3(i,k,j) = cu3(i,k,j)
            end do
        end do
    end do
! Store the scalar in Fourier space
    do n = 1,n_th
        do j = 0,tnky_th
            do k = 0,tnkz_s_th
                do i = 0,nkx_s_th
                    crth(i,k,j,n) = cth(i,k,j,n)
                end do
            end do
        end do
    end do

    if (rank == 0) then
        write(*,*) 'Calling FFT'
    end if

    call mpi_barrier(mpi_comm_world,ierror)

! Now, convert the velocity and vertical gradients to physical space
    call fft_xzy_mpi_to_physical(cu1,u1)
    call fft_xzy_mpi_to_physical(cu2,u2)
    call fft_xzy_mpi_to_physical(cu3,u3)
    do n = 1,n_th
        call fft_xzy_mpi_to_physical_th(cth(0,0,0,n),th(0,0,0,n))
    end do

    if (rank == 0) then
        write(*,*) 'Done FFT'
    end if

    if (movie) then
        ranky_mov = ny_mov/(ny_s+1)
        rankz_mov = nz_mov/(nz_s+1)
        call mpi_barrier(mpi_comm_world,ierror)
        if (ranky==ranky_mov) then
            call writehdf5_xzplane(mod(ny_mov,ny_s+1))
        end if
        call mpi_barrier(mpi_comm_world,ierror)
        if (rankz==rankz_mov) then
            call writehdf5_xyplane(mod(nz_mov,nz_s+1))
        end if
        call mpi_barrier(mpi_comm_world,ierror)
        call writehdf5_yzplane(nx_mov)
    end if

! Get the horizontal means of U1, U3 & TH
    call horiz_mean(cr1,u1me)
    call horiz_mean(cr3,u3me)
    call horiz_mean(crth,th1me)
    call mpi_barrier(mpi_comm_world,ierror)
    if (rank == 0) write(*,*) 'horizontal means computed'

    if (rankz == 0) then
        gname = 'u1me'
        call WriteMeanH5(gname,u1me)
        gname = 'u3me'
        call WriteMeanH5(gname,u3me)
        gname = 'thme'
        call WriteMeanH5(gname,th1me)
    end if
    call mpi_barrier(mpi_comm_world,ierror)

    if (rank == 0) write(*,*) 'horizontal means written to file'
     
! Get the bulk rms values and write perturbation velocities to CRx
    u1rms = 0.d0
    u2rms = 0.d0
    u3rms = 0.d0
    u1u2 = 0.d0
    u3u2 = 0.d0
    u1u1 = 0.d0
    u2u2 = 0.d0
    u3u3 = 0.d0

    do i = 0,nxm
        do j = 0,ny_s
            do k = 0,nz_s
                s1(i,k,j) = u1(i,k,j)-u1me(j)
                u1rms = u1rms+s1(i,k,j)**2.d0
                u1u1(j) = u1u1(j)+s1(i,k,j)**2.d0
                u1u2(j) = u1u2(j)+s1(i,k,j)*u2(i,k,j)
            end do
        end do
    end do
    call fft_xzy_mpi_to_fourier(s1,cr1)
    do i = 0,nxm
        do j = 0,ny_s
            do k = 0,nz_s
                s1(i,k,j) = u2(i,k,j)-dble(cr2(0,0,0))
                u2rms = u2rms+s1(i,k,j)**2.d0
                u2u2(j) = u2u2(j)+s1(i,k,j)**2.d0
                s1(i,k,j) = u3(i,k,j)-u3me(j)
                u3rms = u3rms+s1(i,k,j)**2.d0
                u3u3(j) = u3u3(j)+s1(i,k,j)**2.d0
                u3u2(j) = u3u2(j)+s1(i,k,j)*u2(i,k,j)
            end do
        end do
    end do
    call fft_xzy_mpi_to_fourier(s1,cr3)

    u1rms = u1rms/dble(nx*ny*nz)
    u2rms = u2rms/dble(nx*ny*nz)
    u3rms = u3rms/dble(nx*ny*nz)
    u1u2 = u1u2/dble(nx*nz)
    u3u2 = u3u2/dble(nx*nz)
    u1u1 = u1u1/dble(nx*nz)
    u2u2 = u2u2/dble(nx*nz)
    u3u3 = u3u3/dble(nx*nz)

    call mpi_allreduce(u1rms, u1rms_sum, 1, &
                    mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
    call mpi_allreduce(u2rms, u2rms_sum, 1, &
                    mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
    call mpi_allreduce(u3rms, u3rms_sum, 1, &
                    mpi_double_precision, mpi_sum, mpi_comm_world, ierror)

    call mpi_allreduce(u1u2, u1u2_sum, ny_s+1, &
                    mpi_double_precision, mpi_sum, mpi_comm_z, ierror)
    call mpi_allreduce(u3u2, u3u2_sum, ny_s+1, &
                    mpi_double_precision, mpi_sum, mpi_comm_z, ierror)
    call mpi_allreduce(u1u1, u1rms_h, ny_s+1, &
                    mpi_double_precision, mpi_sum, mpi_comm_z, ierror)
    call mpi_allreduce(u2u2, u2rms_h, ny_s+1, &
                    mpi_double_precision, mpi_sum, mpi_comm_z, ierror)
    call mpi_allreduce(u3u3, u3rms_h, ny_s+1, &
                    mpi_double_precision, mpi_sum, mpi_comm_z, ierror)

    if (rank == 0) then
        u1rms_sum = sqrt(u1rms_sum)
        u2rms_sum = sqrt(u2rms_sum)
        u3rms_sum = sqrt(u3rms_sum)

        write(*,*) '<u1rms>: ',u1rms_sum
        write(*,*) '<u2rms>: ',u2rms_sum
        write(*,*) '<u3rms>: ',u3rms_sum
        
        gname = 'u1rms'
        call WriteStatH5(gname,u1rms_sum)
        gname = 'u2rms'
        call WriteStatH5(gname,u2rms_sum)
        gname = 'u3rms'
        call WriteStatH5(gname,u3rms_sum)
    end if

    if (rankz == 0) then
        gname = 'u1u2'
        call WriteMeanH5(gname,u1u2_sum)
        gname = 'u3u2'
        call WriteMeanH5(gname,u3u2_sum)
        u1rms_h = sqrt(u1rms_h)
        u2rms_h = sqrt(u2rms_h)
        u3rms_h = sqrt(u3rms_h)
        gname = 'u1rms'
        call WriteMeanH5(gname,u1rms_h)
        gname = 'u2rms'
        call WriteMeanH5(gname,u2rms_h)
        gname = 'u3rms'
        call WriteMeanH5(gname,u3rms_h)
    end if
    call mpi_barrier(mpi_comm_world,ierror)

! Do over the number of passive scalars
    do n = 1,n_th

! Get the bacterial/nutrient correlation
        thrms(n) = 0.d0
        thflux(n) = 0.d0
        thu2(:,n) = 0.d0
        do j = 0,ny_s_th
            do k = 0,nz_s_th
                do i = 0,nxm_th
                    sth1(i,k,j) = th(i,k,j,n)-th1me(j)
                    thrms(n) = thrms(n)+sth1(i,k,j)**2
                    thflux(n) = thflux(n)+th(i,k,j,n)*u2(i,k,j)
                    thu2(j,n) = thu2(j,n)+th(i,k,j,n)*u2(i,k,j)
                    thth(j,n) = thth(j,n)+sth1(i,k,j)**2
                end do
            end do
        end do
        call fft_xzy_mpi_to_fourier_th(sth1,crth(:,:,:,1))
        thrms(n) = thrms(n)/dble(nx_th*ny_th*nz_th)
        thflux(n) = thflux(n)/dble(nx_th*ny_th*nz_th)
        thu2(:,n) = thu2(:,n)/dble(nx_th*nz_th)
        thth(:,n) = thth(:,n)/dble(nx_th*nz_th)

        call mpi_allreduce(thrms(n), thrms_sum(n), 1, mpi_double_precision, &
                                mpi_sum, mpi_comm_world, ierror)
        call mpi_allreduce(thflux(n), thflux_sum(n), 1, &
                                mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
        call mpi_allreduce(thu2, thu2_sum, ny_s+1, &
                                mpi_double_precision, mpi_sum, mpi_comm_z, ierror)
        call mpi_allreduce(thth, thth_h, ny_s+1, &
                                mpi_double_precision, mpi_sum, mpi_comm_z, ierror)

        if (rank == 0) then
            thrms_sum = sqrt(thrms_sum)
            thflux_sum = ri_tau(n)*thflux_sum

            write(*,*) 'n,thrms(n),thflux(n): ', &
                        n,thrms_sum(n),thflux_sum(n)

            gname = 'thrms'
            call WriteStatH5(gname,thrms_sum(n))
            gname = 'thflux'
            call WriteStatH5(gname,thflux_sum(n))
        end if

        if (rankz == 0) then
            thu2_sum = Ri_tau(n)*thu2_sum
            gname = 'thflux'
            call WriteMeanH5(gname,thu2_sum(:,n))
            thth_h = sqrt(thth_h)
            gname = 'thrms'
            call WriteMeanH5(gname,thth_h(:,n))
        end if
        call mpi_barrier(mpi_comm_world,ierror)

! End do over number of passive scalars, n
    end do

    call mpi_barrier(mpi_comm_world,ierror)

! Convert back to Fourier space
    do n = 1,n_th
        call fft_xzy_mpi_to_fourier_th(th(0,0,0,n),cth(0,0,0,n))
    end do

    call mpi_barrier(mpi_comm_world,ierror)

! Convert velocity back to Fourier space
    call fft_xzy_mpi_to_fourier(U1,CU1)
    call fft_xzy_mpi_to_fourier(U2,CU2)
    call fft_xzy_mpi_to_fourier(U3,CU3)

end if


call spectra_per(CU1,E)
spectrum = 0.d0
call MPI_ALLREDUCE(E,spectrum,TNKY+1,MPI_DOUBLE_PRECISION, &
                    MPI_SUM,MPI_COMM_WORLD,IERROR)
gname = 'U1'
if (RANK == 0) call WriteSpectrumH5(gname,spectrum)

call spectra_per(CU2,E)
spectrum = 0.d0
call MPI_ALLREDUCE(E,spectrum,TNKY+1,MPI_DOUBLE_PRECISION, &
                    MPI_SUM,MPI_COMM_WORLD,IERROR)
gname = 'U2'
if (RANK == 0) call WriteSpectrumH5(gname,spectrum)

call spectra_per(CU3,E)
spectrum = 0.d0
call MPI_ALLREDUCE(E,spectrum,TNKY+1,MPI_DOUBLE_PRECISION, &
                    MPI_SUM,MPI_COMM_WORLD,IERROR)
gname = 'U3'
if (RANK == 0) call WriteSpectrumH5(gname,spectrum)

call spectra_per(CTH(:,:,:,1),E)
spectrum = 0.d0
call MPI_ALLREDUCE(E,spectrum,TNKY+1,MPI_DOUBLE_PRECISION, &
                    MPI_SUM,MPI_COMM_WORLD,IERROR)
gname = 'TH1'
if (RANK == 0) call WriteSpectrumH5(gname,spectrum)

call tkebudget_per


RETURN
END

!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
subroutine tkebudget_per
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
! Note, it is important to only run this routine after complete R-K
!  time advancement since F1 is overwritten which is needed between R-K steps

include 'header'

real*8 epsilon_mean,epsilon_sum,eta
real*8 chi_mean,chi_sum
real*8 epsilon_h(0:ny_s), chi_h(0:ny_s)
real*8 epsilon_h_sum(0:ny_s), chi_h_sum(0:ny_s)
character(10) :: gname

integer i,j,k

! Compute the scalar dissipation rate, chi = <dth/dx_i dth/dx_i>
chi_mean = 0.d0
chi_h(:) = 0.d0
! Store dth/dx in CSTH1
do j = 0,tnky_th
    do k = 0,tnkz_s_th
        do i = 0,nkx_s_th
            csth1(i,k,j) = cikx_s_th(i)*crth(i,k,j,1)
        end do
    end do
end do
! Convert to physical space
call fft_xzy_mpi_to_physical_th(csth1,sth1)
do j = 0,ny_s_th
    do k = 0,nz_s_th
        do i = 0,nxm_th
            chi_mean = chi_mean+(sth1(i,k,j)**2.0)
            chi_h(j) = chi_h(j)+(sth1(i,k,j)**2.0)
        end do
    end do
end do
! Store dth/dy in CSTH1
do j = 0,tnky_th
    do k = 0,tnkz_s_th
        do i = 0,nkx_s_th
            csth1(i,k,j) = ciky_th(j)*crth(i,k,j,1)
        end do
    end do
end do
! Convert to physical space
call fft_xzy_mpi_to_physical_th(csth1,sth1)
do j = 0,ny_s_th
    do k = 0,nz_s_th
        do i = 0,nxm_th
            chi_mean = chi_mean+(sth1(i,k,j)**2.0)
            chi_h(j) = chi_h(j)+(sth1(i,k,j)**2.0)
        end do
    end do
end do
! Store dth/dz in CSTH1
do j = 0,tnky_th
    do k = 0,tnkz_s_th
        do i = 0,nkx_s_th
            csth1(i,k,j) = cikz_s_th(k)*crth(i,k,j,1)
        end do
    end do
end do
! Convert to physical space
call fft_xzy_mpi_to_physical_th(csth1,sth1)
do j = 0,ny_s_th
    do k = 0,nz_s_th
        do i = 0,nxm_th
            chi_mean = chi_mean+(sth1(i,k,j)**2.0)
            chi_h(j) = chi_h(j)+(sth1(i,k,j)**2.0)
        end do
    end do
end do
! Compute the turbulent dissipation rate, epsilon = nu*<du_i/dx_j du_i/dx_j>
epsilon_mean = 0.d0
epsilon_h = 0.d0
! Store du/dx in CS1
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cs1(i,k,j) = cikx_s(i)*cr1(i,k,j)
        end do
    end do
end do
! Convert to physical space
call fft_xzy_mpi_to_physical(cs1,s1)
do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            epsilon_mean = epsilon_mean+(s1(i,k,j)**2.0)
            epsilon_h(j) = epsilon_h(j)+(s1(i,k,j)**2.0)
        end do
    end do
end do
! Store dv/dx in CS1
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cs1(i,k,j) = cikx_s(i)*cr2(i,k,j)
        end do
    end do
end do
! Convert to physical space
call fft_xzy_mpi_to_physical(cs1,s1)
do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            epsilon_mean = epsilon_mean+s1(i,k,j)**2.0
            epsilon_h(j) = epsilon_h(j)+(s1(i,k,j)**2.0)
        end do
    end do
end do
! Compute du/dy at GYF gridpoints, note remove mean
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cs1(i,k,j) = ciky(j)*cr1(i,k,j)
        end do
    end do
end do
! Convert to physical space
call fft_xzy_mpi_to_physical(cs1,s1)
do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            epsilon_mean = epsilon_mean+s1(i,k,j)**2.0
            epsilon_h(j) = epsilon_h(j)+(s1(i,k,j)**2.0)
        end do
    end do
end do
! Store dw/dx in CS1
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cs1(i,k,j) = cikx_s(i)*cr3(i,k,j)
        end do
    end do
end do
! Convert to physical space
call fft_xzy_mpi_to_physical(cs1,s1)
do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            epsilon_mean = epsilon_mean+s1(i,k,j)**2.0
            epsilon_h(j) = epsilon_h(j)+(s1(i,k,j)**2.0)
        end do
    end do
end do
! Compute du/dz at GYF gridpoints
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cs1(i,k,j) = cikz_s(k)*cr1(i,k,j)
        end do
    end do
end do
! Convert to physical space
call fft_xzy_mpi_to_physical(cs1,s1)
do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            epsilon_mean = epsilon_mean+s1(i,k,j)**2.0
            epsilon_h(j) = epsilon_h(j)+(s1(i,k,j)**2.0)
        end do
    end do
end do
! Compute dv/dy at GYF gridpoints, note remove mean
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cs1(i,k,j) = ciky(j)*cr2(i,k,j)
        end do
    end do
end do
! Convert to physical space
call fft_xzy_mpi_to_physical(cs1,s1)
do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            epsilon_mean = epsilon_mean+s1(i,k,j)**2.0
            epsilon_h(j) = epsilon_h(j)+(s1(i,k,j)**2.0)
        end do
    end do
end do
! Compute dw/dy at GYF gridpoints, note remove mean
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cs1(i,k,j) = ciky(j)*cr3(i,k,j)
        end do
    end do
end do
! Convert to physical space
call fft_xzy_mpi_to_physical(cs1,s1)
do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            epsilon_mean = epsilon_mean+s1(i,k,j)**2.0
            epsilon_h(j) = epsilon_h(j)+(s1(i,k,j)**2.0)
        end do
    end do
end do
! Store dv/dz in CS1
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cs1(i,k,j) = cikz_s(k)*cr2(i,k,j)
        end do
    end do
end do
! Convert to physical space
call fft_xzy_mpi_to_physical(cs1,s1)
do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            epsilon_mean = epsilon_mean+s1(i,k,j)**2.0
            epsilon_h(j) = epsilon_h(j)+(s1(i,k,j)**2.0)
        end do
    end do
end do
! Store dw/dz in CS1
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            cs1(i,k,j) = cikz_s(k)*cr3(i,k,j)
        end do
    end do
end do
! Convert to physical space
call fft_xzy_mpi_to_physical(cs1,s1)
do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            epsilon_mean = epsilon_mean+(s1(i,k,j)**2.0)
            epsilon_h(j) = epsilon_h(j)+(s1(i,k,j)**2.0)
        end do
    end do
end do

if (rank == 0) write(*,*) 'nx,ny,nz: ',nx,ny,nz

epsilon_mean = nu*epsilon_mean/float(nx*ny*nz)
chi_mean = Ri_tau(1)*nu/Pr(1)*chi_mean/float(nx*ny*nz)

call mpi_allreduce(epsilon_mean,epsilon_sum,1, &
                mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
call mpi_allreduce(chi_mean,chi_sum,1, &
                mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

call mpi_allreduce(epsilon_h,epsilon_h_sum,ny_s+1, &
                mpi_double_precision,mpi_sum,mpi_comm_z,ierror)
call mpi_allreduce(chi_h,chi_h_sum,ny_s+1, &
                mpi_double_precision,mpi_sum,mpi_comm_z,ierror)

epsilon_h_sum = nu*epsilon_h_sum/float(nx*nz)
chi_h_sum = Ri_tau(1)*nu/Pr(1)*chi_h_sum/float(nx*nz)

if (rank == 0) write(*,*) 'epsilon_mean: ',epsilon_sum

if (rank == 0) then
    gname = 'epsilon'
    call WriteStatH5(gname,epsilon_sum)
    gname = 'chi'
    call WriteStatH5(gname,chi_sum)
end if

if (rankz == 0) then
    gname = 'epsilon'
    call WriteMeanH5(gname,epsilon_h_sum)
    gname = 'chi'
    call WriteMeanH5(gname,chi_h_sum)
end if
call mpi_barrier(mpi_comm_world,ierror)

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
subroutine spectra_per(cu,E)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
! WORK IN PROGRESS SPECTRUM CALCULATOR
include 'header'

complex*16 cu(0:nx_s/2,0:nz_s,0:ny+1)
real*8 E(0:tnky)
integer i,j,k

E = 0.d0
cs1 = 0.5*cu*conjg(cu)
do j = 0,tnky
    do k = 0,tnkz_s
        do i = 0,nkx_s
            if ((kx_s(i) == 0) .and. ((kz_s(k) < 0) .or. &
                ((kz_s(k) == 0) .and. (ky(j) < 0)))) then
            else
                if ((i+rankz*(nkx_s+1) <= nkx) .and. &
                    (k+ranky*(tnkz_s+1) <= tnkz)) then
                    if ((rank == 0).and. &
                        (i == 0) .and. (j == 0) .and. (k == 0)) then
                        E(j) = E(j)+real(cs1(i,k,j))
                    else
                        E(j) = E(j)+2.d0*real(cs1(i,k,j))
                    end if
                end if
            end if
    
        end do
    end do
end do

end subroutine spectra_per
