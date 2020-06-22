!******************************************************************************|
! fft.f90, the FFT package for diablo.                               VERSION 2.5
!
! This file isolates all calls to the FFTW package (available at: www.fftw.org)
! These wrapper routines were written by T. Bewley (spring 2001).
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
!        do k = 0,tnkz           [where tnkz = 2*nkz = 2*(nz/3) ]
!          do i = 0,nkx          [where  nkx = nx/3             ]
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
subroutine init_fft
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer i,j,k

! Looks like these should be changed to the following form:
! integer, parameter :: fftw_forward = -1
integer :: fftw_forward,        fftw_backward,
integer :: fftw_estimate,       fftw_measure,
integer :: fftw_out_of_place,   fftw_in_place,
integer :: fftw_use_wisdom,     fftw_threadsafe
parameter(  fftw_forward=-1,    fftw_backward=1,&
            fftw_estimate=0,    fftw_measure=1,&
            fftw_out_of_place=0,fftw_in_place=8,&
            fftw_use_wisdom=16, fftw_threadsafe=128)

if (rank.eq.0) then
    write(6,*) 'Initializing FFTW package.'
end if

pi = 4.d0 * atan(1.d0)
ci = cmplx(0.0,1.0)
eps= 0.000000001

call rfftwnd_f77_create_plan(fftw_x_to_f_plan, 1, nx,&
            fftw_forward, fftw_measure)
call rfftwnd_f77_create_plan(fftw_x_to_p_plan, 1, nx,&
            fftw_backward, fftw_measure)
call rfftw_f77_create_plan(fftw_test_plan, nx, fftw_backward, fftw_measure)
nkx = nx/3
rnx = 1.0*nx
do i = 0,nkx
    kx(i) = i*(2.*pi)/lx
    kx2(i) = kx(i)*kx(i)
    cikx(i) = ci*kx(i)
end do

call fftwnd_f77_create_plan(fftw_z_to_f_plan, 1, nz,&
            fftw_forward,  fftw_measure + fftw_in_place)
call fftwnd_f77_create_plan(fftw_z_to_p_plan, 1, nz,&
            fftw_backward, fftw_measure + fftw_in_place)
nkz = nz/3
tnkz = 2*nkz
rnz = 1.0*nz
do k = 0,nkz
    kz(k) = k*(2.*pi)/lz
end do
do k = 1,nkz
    kz(tnkz+1-k) = -k*(2.*pi)/lz
end do
do k = 0,tnkz
    kz2(k) = kz(k)*kz(k)
    cikz(k) = ci*kz(k)
end do

call fftwnd_f77_create_plan(fftw_y_to_f_plan, 1, ny,&
        fftw_forward,  fftw_measure + fftw_in_place)
call fftwnd_f77_create_plan(fftw_y_to_p_plan, 1, ny,&
        fftw_backward, fftw_measure + fftw_in_place)
nky = ny/3
tnky = 2*nky
rny = 1.0*ny
ky(0)  =  0.
ky2(0) = 0.
ciky(0) = (0.0,0.0)
do j = 1,nky
    ky(j) = j*(2.*pi)/ly
end do
do j = 1,nky
    ky(tnky+1-j) = -j*(2.*pi)/ly
end do
do j = 1,tnky
    ky2(j) = ky(j)*ky(j)
    ciky(j) = ci*ky(j)
end do

nkx_s = nint((nkx+1)/np_s+0.5)-1
tnkz_s = nint((tnkz+1)/np_s+0.5)-1

do i = 0,nkx_s
    if ((i+mod(rank,np_s)*(nkx_s+1)).le.nkx) then
        cikx_s(i) = cikx(i+mod(rank,np_s)*(nkx_s+1))
        kx_s(i) = kx(i+mod(rank,np_s)*(nkx_s+1))
        kx2_s(i) = kx2(i+mod(rank,np_s)*(nkx_s+1))
    end if
end do
do k = 0,tnkz_s
    if ((k+int(rank/np_s)*(tnkz_s+1)).le.tnkz) then
        cikz_s(k) = cikz(k+int(rank/np_s)*(tnkz_s+1))
        kz_s(k) = kz(k+int(rank/np_s)*(tnkz_s+1))
        kz2_s(k) = kz2(k+int(rank/np_s)*(tnkz_s+1))
    end if
end do

if (rank.eq.0) then
    write(6,*) 'FFTW package initialized.'
end if

return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_xzy_mpi_to_physical(cu,u)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'

!!! Looks like REAL*8 is deprecated and should use:
! integer, parameter :: dp = kind(1.d0)
! real(dp) :: u(0:nx+1, 0:nz_s, 0:ny_s+1)

integer jmin, jmax, kmin, kmax, i, j, k
! Input array is in Fourier space and local in X
complex*16 cu(0:nx_s/2,0:nz_s,0:ny+1)
! Output array is in physical space and local in Y
real*8     u (0:nx+1,0:nz_s,0:ny_s+1)

! Intermediate arrays local X, and Z, These are all equivalenced
real*8 u_x(0:nx+1,0:nz_s,0:ny_s)
complex*16 cu_x(0:nx/2,0:nz_s,0:ny_s)

! CU_Y and CU_Z can't be equivalenced since they are used
! in the transpose routines
complex*16 cu_y(0:nx_s/2,0:nz_s,0:ny+1)
complex*16 cu_z(0:nx_s/2,0:nz+1,0:ny_s)

! Equivalence the intermediate arrrays to avoid wasting memory
! The FFTs are done in-place, so this is safe
equivalence(cu_x,u_x)

! Inverse transform in the x-direction:
call fft_y_to_physical(cu,cu_y)
call mpi_barrier(mpi_comm_world,ierror)

! Perform a parallel transpose to get data stored locally in the z-direction
call transpose_mpi_y_to_z(cu_y,cu_z)
call mpi_barrier(mpi_comm_world,ierror)

! Inverse transform in the z-direction:
call fft_z_to_physical(cu_z,cu_z)
call mpi_barrier(mpi_comm_world,ierror)

! Perform a parallel transpose to get data stored locally in the y-direction
call transpose_mpi_z_to_x(cu_z,cu_x)
call mpi_barrier(mpi_comm_world,ierror)

! Inverse transform in the y-direction:
call fft_x_to_physical(cu_x,u_x)
call mpi_barrier(mpi_comm_world,ierror)

! Now, transfer from U_X to the output array
do j=0,ny_s
    do k=0,nz_s
        do i=0,nx+1
            u(i,k,j)=u_x(i,k,j)
        end do
    end do
end do

! Set data above NXM=0
do i=nx,nx+1
    u(i,0:nz_s,0:ny_s)=0.d0
end do
! Set data above NYM=0
do j=max(0,ny-int(rank/np_s)*(ny_s+1)),ny_s+1
    u(0:nx+1,0:nz_s,j)=0.d0
end do
! Set data above NZM=0
do k=max(0,nz-mod(rank,np_s)*(nz_s+1)),nz_s
    u(0:nx+1,k,0:ny_s)=0.d0
end do

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_xzy_mpi_to_fourier(u,cu)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms (in 1 direction) planes JMIN-JMAX to Fourier space.
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'

integer jmin, jmax, kmin, kmax, i, j, k
! Input array is in Physical space and local in X
real*8     u (0:nx+1,0:nz_s,0:ny_s+1)
! Output array is in Fourier space and local in Y
complex*16 cu(0:nx_s/2,0:nz_s,0:ny+1)

! Intermediate arrays local X, and Z
real*8 u_x(0:nx+1,0:nz_s,0:ny_s)
complex*16 cu_x(0:nx/2,0:nz_s,0:ny_s)

complex*16 cu_z(0:nkx_s,0:nz+1,0:ny_s)
complex*16 cu_y(0:nkx_s,0:tnkz_s,0:ny+1)

real*8 rnum

equivalence (u_x,cu_x)

do j=0,ny_s
    do k=0,nz_s
        do i=0,nx+1
            u_x(i,k,j)=u(i,k,j)
        end do
    end do
end do

call fft_x_to_fourier(u_x,cu_x)

call mpi_barrier(mpi_comm_world,ierror)

! Parallel transpose to get data stored locally in the z-direction
call transpose_mpi_x_to_z(cu_x,cu_z)
call mpi_barrier(mpi_comm_world,ierror)

! FFT in the z-direction:
call fft_z_to_fourier(cu_z,cu_z)

call mpi_barrier(mpi_comm_world,ierror)

! Parallel transpose to get data stored locally in the z-direction
call transpose_mpi_z_to_y(cu_z,cu_y)
call mpi_barrier(mpi_comm_world,ierror)

! FFT in the y-direction:
call fft_y_to_fourier(cu_y,cu_y)
call mpi_barrier(mpi_comm_world,ierror)

do j=0,ny+1
    do k=0,tnkz_s
        do i=0,nkx_s
            cu(i,k,j)=cu_y(i,k,j)
        end do
        do i=nkx_s+1,nx_s/2
            cu(i,k,j)=0.d0
        end do
    end do
end do

return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_z_to_physical(cu,u)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms along the z direction to Fourier space
! The input and output are real arrays with the input array packed
! to contain both the real and imagingary parts of the transformed array
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer  i, j, k
complex*16   cu(0:nx_s/2,0:nz+1,0:ny_s)
complex*16   u (0:nx_s/2,0:nz+1,0:ny_s)

complex*16 czx_plane(0:nz,0:nx_s/2)

! First, unpack data into the CZX_PLANE complex temporary
! storage variable (dealiasing in the process)
! Then, perform a complex -> complex transform in the z-direction

do j=0,ny_s
    do i=0,nkx_s
        do k=0,nkz
            czx_plane(k,i)=cu(i,k,j)
        end do
        do k=nkz+1,nzm-nkz
            czx_plane(k,i)=cmplx(0.0,0.0)
        end do
        do k=1,nkz
            czx_plane(nzm-nkz+k,i)=cu(i,nkz+k,j)
        end do
    end do
    call fftwnd_f77(fftw_z_to_p_plan, nkx_s+1,&
            czx_plane(0,0), 1, nz+1, czx_plane(0,0), 1, nz+1)

    do k=0,nzm
        do i=0,nkx_s
            u(i,k,j)=czx_plane(k,i)
        end do
        do i=nkx_s+1,nx_s/2
            u(i,k,j)=0.d0
        end do
    end do
end do

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_z_to_fourier(u,cu)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms in the z-direction
! The input and output should be in physical space with the output
! packed to hold the real and imaginary parts
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

include 'header'
integer i, j, k
complex*16 u(0:nkx_s,0:nz+1,0:ny_s)
complex*16 cu(0:nkx_s,0:nz+1,0:ny_s)

complex*16 in_plane(0:nz,0:nkx_s)
complex*16 out_plane(0:nz,0:nkx_s)

! First, put the data into the CZX_PLANE temporary storage variable,
! perform a real -> complex transform in the z direction.


do j=0,ny_s

    do i=0,nkx_s
        do k=0,nzm
            in_plane(k,i)=u(i,k,j)
        end do
    end do


    call fftwnd_f77(fftw_z_to_f_plan, nkx_s+1,&
            in_plane(0,0), 1, nz+1, in_plane(0,0), 1, nz+1)


! Scale by NZ (necessary by FFTW convention)

    do i=0,nkx_s
        do k=0,nkz
            cu(i,k,j)=in_plane(k,i)/rnz
        end do
        do k=1,nkz
            cu(i,nkz+k,j)=in_plane(nzm-nkz+k,i)/rnz
        end do
    end do

end do

return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_y_to_fourier(u,cu)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms in the y-direction
! The input and output should be in physical space with the output
! packed to hold the real and imaginary parts
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer i, j, k
complex*16  u (0:nkx_s,0:tnkz_s,0:ny+1)
complex*16  cu(0:nkx_s,0:tnkz_s,0:ny+1)

complex*16 cyz_plane(0:ny,0:tnkz_s)

! First, put the data into the CYZ_PLANE temporary storage variable,
! perform a real -> complex transform in the y direction.

do i=0,nkx_s
    do k=0,tnkz_s
        do j=0,ny-1
            cyz_plane(j,k)=u(i,k,j)
        end do
    end do

    call fftwnd_f77(fftw_y_to_f_plan, tnkz_s+1,&
                cyz_plane(0,0), 1, ny+1, cyz_plane(0,0), 1, ny+1)

! Scale by NY (necessary by FFTW convention)
    do k=0,tnkz_s
        do j=0,nky
            cu(i,k,j)=cyz_plane(j,k)/rny
        end do
        do j=1,nky
            cu(i,k,nky+j)=cyz_plane(nym-nky+j,k)/rny
        end do
    end do
end do

return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_y_to_physical(cu,u)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms along the y direction to Fourier space
! The input and output are real arrays with the input array packed
! to contain both the real and imagingary parts of the transformed array
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer  i, j, k
complex*16   cu(0:nx_s/2,0:nz_s,0:ny+1)
complex*16   u (0:nx_s/2,0:nz_s,0:ny+1)

complex*16 cyz_plane(0:ny,0:nz_s)

! First, unpack data into the CYZ_PLANE complex temporary
! storage variable (dealiasing in the process)
! Then, perform a complex -> real transform in the y-direction

do i=0,nkx_s
    do k=0,nz_s
! Zero CYZ_PLANE
        do j=0,nky
            cyz_plane(j,k)=cu(i,k,j)
        end do
        do j=nky+1,nym-nky
            cyz_plane(j,k)=cmplx(0.0,0.0)
        end do
        do j=1,nky
            cyz_plane(nym-nky+j,k)=cu(i,k,nky+j)
        end do
    end do
    call fftwnd_f77(fftw_y_to_p_plan, nz_s+1,&
            cyz_plane(0,0), 1, ny+1, cyz_plane(0,0), 1, ny+1)
    do k=0,nz_s
        do j=0,nym
            u(i,k,j)=cyz_plane(j,k)
        end do
    end do
end do

return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_x_to_fourier(u,cu)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms (in 1 direction) to Fourier space.
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
integer i, j, k
real*8     u (0:nx+1,0:nz_s,0:ny_s)
complex*16 cu(0:nx/2,0:nz_s,0:ny_s)

! Looping over the planes of interest, simply perform a real -> complex
! transform in place in the big storage array, scaling appropriately.

do j=0,ny_s
    call rfftwnd_f77_real_to_complex(fftw_x_to_f_plan,(nz_s+1),&
                u(0,0,j), 1, nx+2, cu(0,0,j), 1, nx/2+1)
    do k=0,nz_s
        do i=0,nkx
            cu(i,k,j)=cu(i,k,j)/rnx
        end do
    end do
end do

do j=0,ny_s
    do k=0,nz_s
        do i=nkx+1,nx/2
            cu(i,k,j)=0.d0
        end do
    end do
end do

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine fft_x_to_physical(cu,u)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! This routine transforms (in 1 direction) to physical space.
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'

integer i, j, k
real*8     u (0:nx+1,0:nz_s,0:ny_s)
complex*16 cu(0:nx/2,0:nz_s,0:ny_s)

! Looping over the planes of interest, simply set the higher wavenumbers to
! zero and then perform a complex -> real transform in place in the big
! storage array.

do j=0,ny_s
    do k=0,nz_s
        do i=nkx+1,nx/2
! Perform De-aliasing on the largest third of the wavenumbers
            cu(i,k,j)=0.
        end do
    end do

    call rfftwnd_f77_complex_to_real(fftw_x_to_p_plan,(nz_s+1),&
                            cu(0,0,j), 1, nx/2+1, u(0,0,j), 1, nx+2)

end do


return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine horiz_mean(cu,u_mean)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'

integer j
real*8  u_mean (0:ny_s)
real*8  u_temp (0:ny-1)
complex*16 cu(0:nx_s/2,0:nz_s,0:ny+1)
complex*16 cu_y(0:nx_s/2,0:nz_s,0:ny+1)

if (rank.eq.0) then
    call fft_y_to_physical(cu,cu_y)
    do j=0,ny-1
        u_temp(j)=dble(cu_y(0,0,j))
    end do
end if
call mpi_barrier(mpi_comm_world,ierror)

call mpi_bcast(u_temp, ny, mpi_double_precision, 0, mpi_comm_world, ierror)

do j=0,ny_s
    u_mean(j)=utemp(ranky*(ny_s+1)+j)
end do

end subroutine