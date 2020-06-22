! This file contains subroutines for inputting and outputting data in
! Diablo as well as all subroutines called directly from diablo.f
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine initialize

include 'header'

real    version, current_version
logical reset_time
integer i, j, k, n


open(11, file='input.dat', form='formatted', status='old')

! Read input file.
!   (Note - if you change the following section of code, update the
!    CURRENT_VERSION number to make obsolete previous input files!)

current_version=2.5
read(11,*)
read(11,*)
read(11,*)
read(11,*)
read(11,*) flavor,   version
if (version  /=  current_version) stop 'wrong input data format.'
read(11,*)
read(11,*) nu, lx, ly, lz
read(11,*)
read(11,*) create_new_flow
read(11,*)
read(11,*) n_time_steps, time_limit, delta_t, reset_time, variable_dt, CFL, update_dt
read(11,*)
read(11,*) verbosity, save_flow_int, save_stats_int, movie
read(11,*)
read(11,*) nx_mov, nx_mov_th, ny_mov, ny_mov_th, nz_mov, nz_mov_th
read(11,*)
! Read in the parameters for the N_TH scalars
do n = 1,n_th
    read(11,*)
    read(11,*) create_new_th(n)
    read(11,*)
    read(11,*) Ri_tau(n), Pr(n), reaction(n)
end do

lx=8.d0*atan(1.d0)*lx
ly=8.d0*atan(1.d0)*ly
lz=8.d0*atan(1.d0)*lz

! If we are using MPI, then Initialize the MPI Variables
mpi_io_num=''
mpi_num=0
call init_mpi
if (n_th >= 1) then
    call init_mpi_th
end if

if (movie) then
    call init_movie
    if (rank == 0) write(*,*) 'Movie files created & initialized.'
end if
call init_stats
if (rank == 0) write(*,*) 'Stats file created & initialized.'
call init_mean
if (rank == 0) write(*,*) 'Horiz-mean file created & initialized.'
if (rank == 0) then
    call init_spectra
    call system('mkdir restart_files')
    write(*,*) 'Restart file directory created.'
end if

! Initialize case-specific packages.
call input_per
call create_grid_per
call init_per

! Initialize grid
if ((flavor /= 'Ensemble') .and. (rank == 0)) then
    write(*,*) 'Note that this code is distributed under the GNU General Public License.'
    write(*,*) 'No warranty is expressed or implied.'
    write(*,*)
    write(*,*) 'Flavour: ', flavor
    write(*,*) 'Grid size: nx =', nx,', ny =', ny,', nz =', nz, '.'
    do n = 1,n_th
        write(*,*) 'Scalar number: ', n
        write(*,*) 'Richardson number: ', Ri_tau(n)
        write(*,*) 'Prandtl number: ', Pr(n)
    end do
end if
nxm = nx-1
nym = ny-1
nzm = nz-1
nxm_th = nx_th-1
nym_th = ny_th-1
nzm_th = nz_th-1

! Initialize storage arrays.
do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nx+1
        u1(i,k,j) = 0.
        u3(i,k,j) = 0.
        u2(i,k,j) = 0.
        p (i,k,j) = 0.
        r1(i,k,j) = 0.
        r2(i,k,j) = 0.
        r3(i,k,j) = 0.
        f1(i,k,j) = 0.
        f2(i,k,j) = 0.
        f3(i,k,j) = 0.
        end do
    end do
end do
do j = 0,ny+1
    do k = 0,nz_s
        do i = 0,nx_s/2
        cu1(i,k,j) = 0.
        cu3(i,k,j) = 0.
        cu2(i,k,j) = 0.
        cp (i,k,j) = 0.
        cr1(i,k,j) = 0.
        cr2(i,k,j) = 0.
        cr3(i,k,j) = 0.
        cf1(i,k,j) = 0.
        cf2(i,k,j) = 0.
        cf3(i,k,j) = 0.
        end do
    end do
end do


! Initialize FFT package (includes defining the wavenumber vectors).
write(*,*) 'Initializing FFT', rank
call init_fft
call init_fft_th
call mpi_barrier(mpi_comm_world, ierror)
write(*,*) 'rank, done init fft', rank

! Initialize RKW3 parameters.
h_bar(1) = delta_t*(8.0/15.0)
h_bar(2) = delta_t*(2.0/15.0)
h_bar(3) = delta_t*(5.0/15.0)
beta_bar(1) = 1.0
beta_bar(2) = 25.0/8.0
beta_bar(3) = 9.0/4.0
zeta_bar(1) = 0.0
zeta_bar(2) = -17.0/8.0
zeta_bar(3) = -5.0/4.0

! Initialize values for reading of scalars
num_read_th = 0
do n = 1,n_th
    if (create_new_th(n)) then
        num_read_th = num_read_th
    else
        num_read_th = num_read_th + 1
        read_th_index(num_read_th) = n
    end if
end do
call create_th_per

! Create flow.
if (create_new_flow) then
    call create_flow_per
    if (flavor /= 'Ensemble') then
        write(*,*) 'A new flowfield has been created'
    end if
    if (flavor == 'Basic') then
        write(*,*) 'In diablo_io...'
        call save_stats(.false.)
        call save_flow(.false.)
    end if

else
    call read_flow

! Initialize flow.
    if (reset_time .or. create_new_flow) then
        previous_time_step = 0
        time_step = 0
        time = 0
    end if

    call save_stats(.false.)
    call poisson_p_per

end if

if (flavor == 'Chemotaxis' ) then
    EK = 0.d0
    do j = 0,tnky
        do k = 0,tnkz_s
            do i = 0,nkx_s
                if ((sqrt(kx2_s(i)+ky2(j)+kz2_s(k)) <= 2.5d0) &
                .and.((i+mod(rank,np_s)*(nkx_s+1)) <= nkx) &
                .and.((k+int(rank/np_s)*(tnkz_s+1)) <= tnkz)) then
                    EK = EK + cu1(i,k,j)*conjg(cu1(i,k,j)) + cu2(i,k,j)*conjg(cu2(i,k,j)) &
                        + cu3(i,k,j)*conjg(cu3(i,k,j))
                end if
            end do
        end do
    end do
    call mpi_allreduce(EK, EK_sum, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
    EK0 = EK_sum
    if (rank == 0) write(*,*) 'Initial EK0: ',EK0
end if

return
end



!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine save_stats(final)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
logical final

call save_stats_per(final)

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine read_flow
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
character*60 fname
integer i, j, k, n

fname = 'start.h5'
if (rank.eq.0) write(*,*) 'Reading flow from ', fname

call mpi_barrier(mpi_comm_world, ierror)
call readhdf5(fname)

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine save_flow(final)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'
character*60 fname
integer      i, j, k, n
logical      final, save_pressure

save_pressure = .false.
if (final) then
    fname = 'end.h5'
    save_pressure = .true.
else
    if (save_flow_int > 0) then
        if (n_time_steps/save_flow_int < 100) then
            fname = 'out'//&
                    char(floor(time_step/save_flow_int/10.)+48)//&
                    char(mod(time_step/int(save_flow_int),10)+48)//'.h5'
        else if (n_time_steps/save_flow_int < 1000) then
            fname = 'out'//char(floor(time_step/save_flow_int/100.)+48)//&
                    char(floor(mod(time_step/int(save_flow_int),100)/10.)+48)//&
                    char(mod(time_step/int(save_flow_int),10)+48)//'.h5'
        end if
    else
        fname = 'out'//&
                char(floor(-(time+delta_t)/save_flow_int/10.)+48)//&
                char(mod(floor(-(time+delta_t)/save_flow_int),10)+48)//'.h5'
    end if
end if

call mpi_barrier(mpi_comm_world, ierror)
call writeHDF5('restart_files/'//fname, save_pressure)

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine end_run(flag)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
include 'header'

logical flag

flag = .false.
! Check for the time
call wall_time(end_time)
if (end_time - start_time > time_limit) then
    if (rank.eq.0) then
        write(*,*) 'STOP because of wall-time hit!'
    end if
    flag = .true.
end if

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
subroutine wall_time(wt)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!
! Return wall-clock time as seconds after Jan. 1, 2018.
! Support for leap year is not included anymore.
!
! By using a 'save' statement, the wall-time after the first
! call to the subroutine could be computed, but that is not
! intended with the present subroutine (e.g. the history file)
!
implicit none

real*8 wt
integer val(8), i, shift, day

integer mon(12,2)
data mon /31,28,31,30,31,30,31,31,30,31,30,31,31,29,31,30,31,30,31,31,30,31,30,31/
!
! Get current date and time
! val(1) : year
! val(2) : month
! val(3) : day
! val(4) : difference to GMT
! val(5) : hour
! val(6) : minute
! val(7) : second
! val(8) : 1/1000 second
!
call date_and_time(values = val)
!
! Determine leap year
!
if (mod(val(1),4) == 0) then
    if (mod(val(1),100) == 0) then
        if (mod(val(1),400) == 0) then
            shift = 2
        else
            shift = 1
        end if
    else
        shift = 2
    end if
else
    shift = 1
end if
!
! Construct day of the year
!
day = val(3)-1
do i = 1,val(2)-1
    day = day+mon(i,shift)
end do
!
! And compute wall-clock time
!
wt = (val(1)-2018)*365*86400 + &
    day*86400 + val(5)*3600 + val(6)*60 + val(7) + dble(val(8)/1000.d0)

end