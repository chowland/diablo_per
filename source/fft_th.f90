!******************************************************************************|
! fft_th.f90, the FFT package for diablo.                               VERSION 0.9
!
! This file isolates all calls to the FFTW package (available at: www.fftw.org)
! These wrapper routines were written by T. Bewley (spring 2001).
! This subroutine has been modified to operate on the scalar grid only
!******************************************************************************|

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! The arrangement of the significant real numbers in the arrays (denoted by +)
! in physical space, in Fourier space, and in Fourier space after packing are
! shown below for the 2D (X-Z) plane.  The third direction (Y) is handled in
! an identical matter as the Z direction shown here.
!
!      oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
!      oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
! nz-1 ++++++++++++++++oo     -1  ++++++++++++oooooo         oooooooooooooooooo
!      ++++++++++++++++oo     -2  ++++++++++++oooooo         oooooooooooooooooo
!      ++++++++++++++++oo     -3  ++++++++++++oooooo         oooooooooooooooooo
!      ++++++++++++++++oo         ++++++++++++oooooo         oooooooooooooooooo
!      ++++++++++++++++oo    -nkz ++++++++++++oooooo         oooooooooooooooooo
!      ++++++++++++++++oo         oooooooooooooooooo     -1  ++++++++++++oooooo
!      ++++++++++++++++oo         oooooooooooooooooo     -2  ++++++++++++oooooo
!      ++++++++++++++++oo         oooooooooooooooooo     -3  ++++++++++++oooooo
!      ++++++++++++++++oo         oooooooooooooooooo         ++++++++++++oooooo
!      ++++++++++++++++oo         oooooooooooooooooo    -nkz ++++++++++++oooooo
!      ++++++++++++++++oo     nkz ++++++++++++oooooo     nkz ++++++++++++oooooo
!      ++++++++++++++++oo         ++++++++++++oooooo         ++++++++++++oooooo
!   3  ++++++++++++++++oo      3  ++++++++++++oooooo      3  ++++++++++++oooooo
!   2  ++++++++++++++++oo      2  ++++++++++++oooooo      2  ++++++++++++oooooo
!   1  ++++++++++++++++oo      1  ++++++++++++oooooo      1  ++++++++++++oooooo
!   0  ++++++++++++++++oo      0  +o++++++++++oooooo      0  +o++++++++++oooooo
!      ^^^^           ^           ^ ^ ^     ^                ^ ^ ^     ^
!      0123           nx-1        0 1 2     nkx              0 1 2     nkx
!
!       PHYSICAL SPACE              FOURIER SPACE         FOURIER SPACE (PACKED)
!
! After the Real->Fourier transform, the significant coefficients are put next
! to each other in the array, so a loop such as
!
!        do k=0,tnkz           [where tnkz = 2*nkz = 2*(nz/3) ]
!          do i=0,nkx          [where  nkx = nx/3             ]
!            cp(i,k,j)= ...
!          end do
!        end do
!
! includes all the Fourier coefficients of interest.  The subsequent loops in
! Fourier space just work on these coefficients in the matrix.
!
! Before a Fourier->Real transform, the significant coefficients are unpacked
! and the higher wavenumbers are SET TO ZERO before the inverse transform.
! This has the effect of doing the required dealiasing.
!
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine init_fft_th
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

include 'header'
integer i,j,k

integer :: fftw_forward,        fftw_backward,
integer :: fftw_estimate,       fftw_measure,
integer :: fftw_out_of_place,   fftw_in_place,
integer :: fftw_use_wisdom,     fftw_threadsafe
parameter(  fftw_forward=-1,    fftw_backward=1,&
            fftw_estimate=0,    fftw_measure=1,&
            fftw_out_of_place=0,fftw_in_place=8,&
            fftw_use_wisdom=16, fftw_threadsafe=128)

if (rank == 0) then
    write(6,*) 'Initializing FFTW package for scalar.'
end if

pi = 4.d0 * atan(1.d0)
ci = cmplx(0.0,1.0)
eps= 0.000000001

call rfftwnd_f77_create_plan(fftw_x_to_f_plan_th, 1, nx_th,&
            fftw_forward,  fftw_measure )
call rfftwnd_f77_create_plan(fftw_x_to_p_plan_th, 1, nx_th,&
            fftw_backward, fftw_measure )
call rfftw_f77_create_plan(fftw_test_plan_th,nx_th,fftw_backward,fftw_measure )
nkx_th = nx_th/3
rnx_th = 1.0*nx_th
do i = 0,nkx_th
    kx_th(i) = i*(2.*pi)/lx
    kx2_th(i) = kx_th(i)*kx_th(i)
    cikx_th(i) = ci*kx_th(i)
end do

call fftwnd_f77_create_plan(fftw_z_to_f_plan_th, 1, nz_th,&
            fftw_forward,  fftw_measure + fftw_in_place)
call fftwnd_f77_create_plan(fftw_z_to_p_plan_th, 1, nz_th,&
            fftw_backward, fftw_measure + fftw_in_place)
nkz_th = nz_th/3
tnkz_th = 2*nkz_th
rnz_th = 1.0*nz_th
do k = 0,nkz_th
    kz_th(k) = k*(2.*pi)/lz
end do
do k = 1,nkz_th
    kz_th(tnkz_th+1-k) = -k*(2.*pi)/lz
end do
do k = 0,tnkz_th
    kz2_th(k) = kz_th(k)*kz_th(k)
    cikz_th(k) = ci*kz_th(k)
end do

call fftwnd_f77_create_plan(fftw_y_to_f_plan_th, 1, ny_th,&
            fftw_forward,  fftw_measure + fftw_in_place)
call fftwnd_f77_create_plan(fftw_y_to_p_plan_th, 1, ny_th,&
            fftw_backward, fftw_measure + fftw_in_place)
nky_th = ny_th/3
tnky_th = 2*nky_th
rny_th = 1.0*ny_th
ky_th(0)  =  0.
ky2_th(0)  =  0.
ciky_th(0)  =  (0.0,0.0)
do j = 1,nky_th
    ky_th(j) = j*(2.*pi)/ly
end do
do j = 1,nky_th
    ky_th(tnky_th+1-j) = -j*(2.*pi)/ly
end do
do j = 1,tnky_th
    ky2_th(j) = ky_th(j)*ky_th(j)
    ciky_th(j) = ci*ky_th(j)
end do

nkx_s_th = nint((nkx_th+1)/np_s+0.5)-1
tnkz_s_th = nint((tnkz_th+1)/np_s+0.5)-1

do i = 0,nkx_s_th
    if ((i+mod(rank,np_s)*(nkx_s_th+1)) <= nkx_th) then
        cikx_s_th(i) = cikx_th(i+mod(rank,np_s)*(nkx_s_th+1))
        kx_s_th(i) = kx_th(i+mod(rank,np_s)*(nkx_s_th+1))
        kx2_s_th(i) = kx2_th(i+mod(rank,np_s)*(nkx_s_th+1))
    end if
end do
do k = 0,tnkz_s_th
    if ((k+int(rank/np_s)*(tnkz_s_th+1)) <= tnkz_th) then
        cikz_s_th(k) = cikz_th(k+int(rank/np_s)*(tnkz_s_th+1))
        kz_s_th(k) = kz_th(k+int(rank/np_s)*(tnkz_s_th+1))
        kz2_s_th(k) = kz2_th(k+int(rank/np_s)*(tnkz_s_th+1))
    end if
end do

if (rank == 0) then
    write(6,*) 'FFTW package initialized for scalar.'
end if

return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_xzy_mpi_to_physical_th(cu,u)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'

integer jmin, jmax, kmin, kmax, i, j, k
! Input array is in Fourier space and local in X
complex*16 cu(0:nx_s_th/2,0:nz_s_th,0:ny_th+1)
! Output array is in physical space and local in Y
real*8     u (0:nx_th+1,0:nz_s_th,0:ny_s_th+1)

! Intermediate arrays local X, and Z, These are all equivalenced
real*8 u_x(0:nx_th+1,0:nz_s_th,0:ny_s_th)
complex*16 cu_x(0:nx_th/2,0:nz_s_th,0:ny_s_th)

! CU_Y and CU_Z can't be equivalenced since they are used
! in the transpose routines
complex*16 cu_y(0:nx_s_th/2,0:nz_s_th,0:ny_th+1)
complex*16 cu_z(0:nx_s_th/2,0:nz_th+1,0:ny_s_th)

! Equivalence the intermediate arrrays to avoid wasting memory
! The FFTs are done in-place, so this is safe
equivalence(cu_x,u_x)

! Inverse transform in the x-direction:
call fft_y_to_physical_th(cu,cu_y)
call mpi_barrier(mpi_comm_world,ierror)

! Perform a parallel transpose to get data stored locally in the z-direction
call transpose_mpi_y_to_z_th(cu_y,cu_z)
call mpi_barrier(mpi_comm_world,ierror)

! Inverse transform in the z-direction:
call fft_z_to_physical_th(cu_z,cu_z)
call mpi_barrier(mpi_comm_world,ierror)

! Perform a parallel transpose to get data stored locally in the y-direction
call transpose_mpi_z_to_x_th(cu_z,cu_x)
call mpi_barrier(mpi_comm_world,ierror)

! Inverse transform in the y-direction:
call fft_x_to_physical_th(cu_x,u_x)
call mpi_barrier(mpi_comm_world,ierror)

! Now, transfer from U_X to the output array
do j = 0,ny_s_th
    do k = 0,nz_s_th
        do i = 0,nx_th+1
            u(i,k,j) = u_x(i,k,j)
        end do
    end do
end do

! Set data above NXM = 0
do i = nx_th,nx_th+1
    u(i,0:nz_s,0:ny_s) = 0.d0
end do
! Set data above NYM = 0
do j = max(0,ny_th-int(rank/np_s)*(ny_s_th+1)),ny_s_th+1
    u(0:nx_th+1,0:nz_s_th,j) = 0.d0
end do
! Set data above NZM = 0
do k = max(0,nz_th-mod(rank,np_s)*(nz_s_th+1)),nz_s_th
    u(0:nx_th+1,k,0:ny_s_th) = 0.d0
end do

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_xzy_mpi_to_fourier_th(u,cu)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'

integer jmin, jmax, kmin, kmax, i, j, k
! Input array is in Physical space and local in X
real*8     u (0:nx_th+1,0:nz_s_th,0:ny_s_th+1)
! Output array is in Fourier space and local in Y
complex*16 cu(0:nx_s_th/2,0:nz_s_th,0:ny_th+1)

! Intermediate arrays local X, and Z
real*8 u_x(0:nx_th+1,0:nz_s_th,0:ny_s_th)
complex*16 cu_x(0:nx_th/2,0:nz_s_th,0:ny_s_th)

complex*16 cu_z(0:nkx_s_th,0:nz_th+1,0:ny_s_th)
complex*16 cu_y(0:nkx_s_th,0:tnkz_s_th,0:ny_th+1)

real*8 rnum

equivalence (u_x,cu_x)

do j = 0,ny_s_th
    do k = 0,nz_s_th
        do i = 0,nx_th+1
            u_x(i,k,j) = u(i,k,j)
        end do
    end do
end do

call fft_x_to_fourier_th(u_x,cu_x)

call mpi_barrier(mpi_comm_world,ierror)

! Parallel transpose to get data stored locally in the z-direction
call transpose_mpi_x_to_z_th(cu_x,cu_z)
call mpi_barrier(mpi_comm_world,ierror)

! FFT in the z-direction:
call fft_z_to_fourier_th(cu_z,cu_z)

call mpi_barrier(mpi_comm_world,ierror)

! Parallel transpose to get data stored locally in the z-direction
call transpose_mpi_z_to_y_th(cu_z,cu_y)
call mpi_barrier(mpi_comm_world,ierror)

! FFT in the y-direction:
call fft_y_to_fourier_th(cu_y,cu_y)
call mpi_barrier(mpi_comm_world,ierror)

do j = 0,ny_th+1
    do k = 0,tnkz_s_th
        do i = 0,nkx_s_th
            cu(i,k,j) = cu_y(i,k,j)
        end do
        do i = nkx_s_th+1,nx_s_th/2
            cu(i,k,j) = 0.d0
        end do
    end do
end do

return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_z_to_physical_th(cu,u)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms along the z direction to Fourier space
! The input and output are real arrays with the input array packed
! to contain both the real and imagingary parts of the transformed array
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer  i, j, k
complex*16   cu(0:nx_s_th/2,0:nz_th+1,0:ny_s_th)
complex*16   u (0:nx_s_th/2,0:nz_th+1,0:ny_s_th)

complex*16 czx_plane(0:nz_th,0:nx_s_th/2)

! First, unpack data into the CZX_PLANE complex temporary
! storage variable (dealiasing in the process)
! Then, perform a complex -> complex transform in the z-direction

do j = 0,ny_s_th
    do i = 0,nkx_s_th
        do k = 0,nkz_th
            czx_plane(k,i) = cu(i,k,j)
        end do
        do k = nkz_th+1,nzm_th-nkz_th
            czx_plane(k,i) = cmplx(0.0,0.0)
        end do
        do k = 1,nkz_th
            czx_plane(nzm_th-nkz_th+k,i) = cu(i,nkz_th+k,j)
        end do
    end do
    call fftwnd_f77(fftw_z_to_p_plan_th, nkx_s_th+1,&
                czx_plane(0,0), 1, nz_th+1, czx_plane(0,0), 1, nz_th+1)

    do k = 0,nzm_th
        do i = 0,nkx_s_th
            u(i,k,j) = czx_plane(k,i)
        end do
        do i = nkx_s_th+1,nx_s_th/2
            u(i,k,j) = 0.d0
        end do
    end do
end do

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_z_to_fourier_th(u,cu)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms in the z-direction
! The input and output should be in physical space with the output
! packed to hold the real and imaginary parts
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

include 'header'
integer i, j, k
complex*16 u(0:nkx_s_th,0:nz_th+1,0:ny_s_th)
complex*16 cu(0:nkx_s_th,0:nz_th+1,0:ny_s_th)

complex*16 in_plane(0:nz_th,0:nkx_s_th)
complex*16 out_plane(0:nz_th,0:nkx_s_th)

! First, put the data into the CZX_PLANE temporary storage variable,
! perform a real -> complex transform in the z direction.


do j = 0,ny_s_th

    do i = 0,nkx_s_th
        do k = 0,nzm_th
            in_plane(k,i) = u(i,k,j)
        end do
    end do


    call fftwnd_f77(fftw_z_to_f_plan_th, nkx_s_th+1,&
                in_plane(0,0), 1, nz_th+1, in_plane(0,0), 1, nz_th+1)


! Scale by NZ (necessary by FFTW convention)

    do i = 0,nkx_s_th
        do k = 0,nkz_th
            cu(i,k,j) = in_plane(k,i)/rnz_th
        end do
        do k = 1,nkz_th
            cu(i,nkz_th+k,j) = in_plane(nzm_th-nkz_th+k,i)/rnz_th
        end do
    end do

end do

return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_y_to_fourier_th(u,cu)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms in the y-direction
! The input and output should be in physical space with the output
! packed to hold the real and imaginary parts
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer i, j, k
complex*16  u (0:nkx_s_th,0:tnkz_s_th,0:ny_th+1)
complex*16  cu(0:nkx_s_th,0:tnkz_s_th,0:ny_th+1)

complex*16 cyz_plane(0:ny_th,0:tnkz_s_th)

! First, put the data into the CYZ_PLANE temporary storage variable,
! perform a real -> complex transform in the y direction.

do i = 0,nkx_s_th
    do k = 0,tnkz_s_th
        do j = 0,ny_th-1
            cyz_plane(j,k) = u(i,k,j)
        end do
    end do

    call fftwnd_f77(fftw_y_to_f_plan_th, tnkz_s_th+1,&
                cyz_plane(0,0), 1, ny_th+1, cyz_plane(0,0), 1, ny_th+1)

! Scale by NY (necessary by FFTW convention)
    do k = 0,tnkz_s_th
        do j = 0,nky_th
            cu(i,k,j) = cyz_plane(j,k)/rny_th
        end do
        do j = 1,nky_th
            cu(i,k,nky_th+j) = cyz_plane(nym_th-nky_th+j,k)/rny_th
        end do
    end do
end do

return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_y_to_physical_th(cu,u)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms along the y direction to Fourier space
! The input and output are real arrays with the input array packed
! to contain both the real and imagingary parts of the transformed array
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer  i, j, k
complex*16   cu(0:nx_s_th/2,0:nz_s_th,0:ny_th+1)
complex*16   u (0:nx_s_th/2,0:nz_s_th,0:ny_th+1)

complex*16 cyz_plane(0:ny_th,0:nz_s_th)

! First, unpack data into the CYZ_PLANE complex temporary
! storage variable (dealiasing in the process)
! Then, perform a complex -> real transform in the y-direction

cyz_plane(:,:) = (0.d0,0.d0)

do i = 0,nkx_s_th
    do k = 0,nz_s_th
        do j = 0,nky_th
            cyz_plane(j,k) = cu(i,k,j)
        end do
        do j = nky_th+1,nym_th-nky_th
            cyz_plane(j,k) = cmplx(0.0,0.0)
        end do
        do j = 1,nky_th
            cyz_plane(nym_th-nky_th+j,k) = cu(i,k,nky_th+j)
        end do
    end do
    call fftwnd_f77(fftw_y_to_p_plan_th, nz_s_th+1,&
                cyz_plane(0,0), 1, ny_th+1, cyz_plane(0,0), 1, ny_th+1)
    do k = 0,nz_s_th
        do j = 0,nym_th
            u(i,k,j) = cyz_plane(j,k)
        end do
    end do
end do

return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_x_to_fourier_th(u,cu)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms (in 1 direction) to Fourier space.
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer i, j, k
real*8     u (0:nx_th+1,0:nz_s_th,0:ny_s_th)
complex*16 cu(0:nx_th/2,0:nz_s_th,0:ny_s_th)

! Looping over the planes of interest, simply perform a real -> complex
! transform in place in the big storage array, scaling appropriately.

do j = 0,ny_s_th
    call rfftwnd_f77_real_to_complex(fftw_x_to_f_plan_th,(nz_s_th+1),&
                    u(0,0,j), 1, nx_th+2, cu(0,0,j), 1, nx_th/2+1)
    do k = 0,nz_s_th
        do i = 0,nkx_th
            cu(i,k,j) = cu(i,k,j)/rnx_th
        end do
    end do
end do

do j = 0,ny_s_th
    do k = 0,nz_s_th
        do i = nkx_th+1,nx_th/2
            cu(i,k,j) = 0.d0
        end do
    end do
end do

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_x_to_physical_th(cu,u)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms (in 1 direction) to physical space.
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'

integer i, j, k
real*8     u (0:nx_th+1,0:nz_s_th,0:ny_s_th)
complex*16 cu(0:nx_th/2,0:nz_s_th,0:ny_s_th)

! Looping over the planes of interest, simply set the higher wavenumbers to
! zero and then perform a complex -> real transform in place in the big
! storage array.

do j = 0,ny_s_th
    do k = 0,nz_s_th
        do i = nkx_th+1,nx_th/2
! Perform De-aliasing on the largest third of the wavenumbers
            cu(i,k,j) = 0.
        end do
    end do

    call rfftwnd_f77_complex_to_real(fftw_x_to_p_plan_th,(nz_s_th+1),&
                cu(0,0,j), 1, nx_th/2+1, u(0,0,j), 1, nx_th+2)

end do


return
end


! The following subroutine is for taking a velocity field in Fourier space
! and interpolating it to the scalar (TH) grid in Physical space
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_xzy_mpi_to_physical_interp(cu,u)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'

integer jmin, jmax, kmin, kmax, i, j, k
! Input array is in Fourier space and local in X on the velocity grid
complex*16 cu(0:nx_s/2,0:nz_s,0:ny+1)
! Output array is in physical space and local in Y on the TH grid
real*8  u (0:nx_th+1,0:nz_s_th,0:ny_s_th+1)

! Intermediate arrays local X, and Z, These are all equivalenced
complex*16 cu_x(0:nx/2,0:nz_s_th,0:ny_s_th)
complex*16 cu_x_th(0:nx_th/2,0:nz_s_th,0:ny_s_th)
real*8 u_x(0:nx_th+1,0:nz_s_th,0:ny_s_th)

! CU_Y and CU_Z can't be equivalenced since they are used
! in the transpose routines
complex*16 cu_y(0:nx_s/2,0:nz_s,0:ny_th+1)

complex*16 cu_z(0:nx_s/2,0:nz+1,0:ny_s_th)

complex*16 cu_z_th(0:nx_s/2,0:nz_th+1,0:ny_s_th)

! Equivalence the intermediate arrrays to avoid wasting memory
! The FFTs are done in-place, so this is safe
equivalence(cu_x_th,u_x)


! First, pad the velocity array with zeros
cu_y(:,:,:) = 0.d0
do j = 0,nky
    cu_y(:,:,j) = cu(:,:,j)
end do
do j = 0,nky-1
    cu_y(:,:,tnky_th-j) = cu(:,:,tnky-j)
end do
! Inverse transform in the x-direction:
call fft_y_to_physical_interp(cu_y,cu_y)
call mpi_barrier(mpi_comm_world,ierror)

! Perform a parallel transpose to get data stored locally in the z-direction
call transpose_mpi_y_to_z_interp(cu_y,cu_z)
call mpi_barrier(mpi_comm_world,ierror)

! Pad the high resolution array with zeros
cu_z_th(:,:,:) = 0.d0
do k = 0,nkz
    cu_z_th(:,k,:) = cu_z(:,k,:)
end do
do k = 0,nkz-1
    cu_z_th(:,tnkz_th-k,:) = cu_z(:,tnkz-k,:)
end do

! Inverse transform in the z-direction:
call fft_z_to_physical_interp(cu_z_th,cu_z_th)
call mpi_barrier(mpi_comm_world,ierror)

! Perform a parallel transpose to get data stored locally in the y-direction
call transpose_mpi_z_to_x_interp(cu_z_th,cu_x)
call mpi_barrier(mpi_comm_world,ierror)

! Pad the high resolution array with zeros
cu_x_th(:,:,:) = 0.d0
do i = 0,nkx
    cu_x_th(i,:,:) = cu_x(i,:,:)
end do

! Inverse transform in the y-direction:
call fft_x_to_physical_th(cu_x_th,u_x)
call mpi_barrier(mpi_comm_world,ierror)

! Now, transfer from U_X to the output array
do j = 0,ny_s_th
    do k = 0,nz_s_th
        do i = 0,nx_th+1
            u(i,k,j) = u_x(i,k,j)
        end do
    end do
end do

! Set data above NXM = 0
do i = nx_th,nx_th+1
    u(i,0:nz_s,0:ny_s) = 0.d0
end do
! Set data above NYM = 0
do j = ny_th-int(rank/np_s)*(ny_s_th+1),ny_s_th+1
    u(0:nx_th+1,0:nz_s_th,j) = 0.d0
end do
! Set data above NZM = 0
do k = nz_th-mod(rank,np_s)*(nz_s_th+1),nz_s_th
    u(0:nx_th+1,k,0:ny_s_th) = 0.d0
end do

return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_z_to_physical_interp(cu,u)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms along the z direction to Fourier space
! The input and output are real arrays with the input array packed
! to contain both the real and imagingary parts of the transformed array
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer  i, j, k
complex*16   cu(0:nx_s/2,0:nz_th+1,0:ny_s_th)
complex*16   u (0:nx_s/2,0:nz_th+1,0:ny_s_th)

complex*16 czx_plane(0:nz_th,0:nx_s/2)

! First, unpack data into the CZX_PLANE complex temporary
! storage variable (dealiasing in the process)
! Then, perform a complex -> complex transform in the z-direction

do j = 0,ny_s_th
    do i = 0,nkx_s
        do k = 0,nkz_th
            czx_plane(k,i) = cu(i,k,j)
        end do
        do k = nkz_th+1,nzm_th-nkz_th
            czx_plane(k,i) = cmplx(0.0,0.0)
        end do
        do k = 1,nkz_th
            czx_plane(nzm_th-nkz_th+k,i) = cu(i,nkz_th+k,j)
        end do
    end do
    call fftwnd_f77(fftw_z_to_p_plan_th, nkx_s+1,&
            czx_plane(0,0), 1, nz_th+1, czx_plane(0,0), 1, nz_th+1)
    do k = 0,nzm_th
        do i = 0,nkx_s
            u(i,k,j) = czx_plane(k,i)
        end do
        do i = nkx_s+1,nx_s/2
            u(i,k,j) = 0.d0
        end do
    end do
end do

return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_y_to_physical_interp(cu,u)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms along the y direction to Fourier space
! The input and output are real arrays with the input array packed
! to contain both the real and imagingary parts of the transformed array
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer  i, j, k
complex*16   cu(0:nx_s/2,0:nz_s,0:ny_th+1)
complex*16   u (0:nx_s/2,0:nz_s,0:ny_th+1)

complex*16 cyz_plane(0:ny_th,0:nz_s)

! First, unpack data into the CYZ_PLANE complex temporary
! storage variable (dealiasing in the process)
! Then, perform a complex -> real transform in the y-direction

cyz_plane(:,:) = (0.d0,0.d0)

do i = 0,nkx_s
    do k = 0,nz_s
        do j = 0,nky_th
            cyz_plane(j,k) = cu(i,k,j)
        end do
        do j = nky_th+1,nym_th-nky_th
            cyz_plane(j,k) = cmplx(0.0,0.0)
        end do
        do j = 1,nky_th
            cyz_plane(nym_th-nky_th+j,k) = cu(i,k,nky_th+j)
        end do
    end do
    call fftwnd_f77(fftw_y_to_p_plan_th, nz_s+1,&
            cyz_plane(0,0), 1, ny_th+1, cyz_plane(0,0), 1, ny_th+1)
    do k = 0,nz_s
        do j = 0,nym_th
            u(i,k,j) = cyz_plane(j,k)
        end do
    end do
end do

return
end