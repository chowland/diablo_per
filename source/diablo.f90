!******************************************************************************|
! diablo.f90 -> DNS In A Box, Laptop Optimized                       VERSION 2.5
!
! This Fortran 90 code computes incompressible flow in a box.
!
! Primative variables (u,v,w,p) are used, and continuity is enforced with a
! fractional step algorithm.
!
! SPATIAL DERIVATIVES:
!   0, 1, 2, or 3 directions are taken to be periodic and handled spectrally
!   (these cases are referred to as the "periodic", "channel", "duct", and
!    "cavity" cases respectively).
!   The remaining directions are taken to be bounded by walls and handled with
!   momentum- and energy-conserving second-order central finite differences.
!
! TIME ADVANCEMENT
!   Two main approaches are implemented:
!     1. RKW3 on nonlinear terms and CN on viscous terms over each RK substep.
!     2. RKW3 on y-derivative terms and CN on other terms over each RK substep.
!
! The emphasis in this introductory code is on code simplicity:
!   -> All variables are in core.
!   -> The code is not explicitly designed for use with either MPI or SMP.
!   -> Overindexing is not used.
! A few simple high-performance programming constructs are used:
!   -> The inner 2 loops are broken out in such a way as to enable out-of-order
!      execution of these loops as much as possible, thereby leveraging
!      vector and superscalar CPU architectures.
!   -> The outer loops are fairly long (including as many operations as
!      possible inside on a single J plane of data) in order to make effective
!      use of cache.
! Multiple time advancement algorithms are implemented for the periodic,
! channel, duct, and cavity cases in order to compare their efficiency for
! various flows on different computational architectures.  In a few of the
! algorithms, the code is broken out fully into a PASS1/PASS2 architecture
! to maximize the efficient use of cache.
!
! This code was developed as a joint project in MAE 223 (CFD), taught by
! Thomas Bewley, at UC San Diego (spring 2001, spring 2005).
! Primary contributions follow:
! Thomas Bewley was the chief software architect
! John R. Taylor wrote the channel flow solvers
!*****************************************************************************|
!
! This code is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation; either version 2 of the License, or (at your
! option) any later version. This code is distributed in the hope that it
! will be useful, but WITHOUT ANY WARRANTY; without even the implied
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details. You should have received a
! copy of the GNU General Public License along with this code; if not,
! write to the Free Software Foundation, Inc., 59 Temple Place - Suite
! 330, Boston, MA 02111-1307, USA.
!
!******************************************************************************|

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
program DIABLO

include 'header'

integer n
logical flag

call initialize

! Initialize START_TIME for run timing
call wall_time(start_time)

if (rank == 0) then
    write(*,*)
    write(*,*) '             ****** Welcome to DIABLO ******'
    write(*,*)
    write(*,*) 'MPI initialized with ', nprocs, ' processes'
end if

! A flag to determine if we are considering the first time-step
first_time = .true.

do time_step = time_step + 1, time_step + n_time_steps
    if (rank == 0) write(*,*) 'now beginning time_step = ',time_step

    do rk_step = 1,3
        call rk_per_1
    end do

    time = time + delta_t
    first_time = .false.

! Save statistics to an output file
    if (save_stats_int < 0) then
        if ((mod(time, abs(save_stats_int)) + delta_t) > abs(save_stats_int)) then
            call save_stats(.false.)
        end if
    else
        if (mod(time_step, int(save_stats_int)) == 0) then
            call save_stats(.false.)
        end if
    end if
! Save the flow to a restart file
    if (save_flow_int < 0) then
        if ((mod(time, abs(save_flow_int)) + delta_t) > abs(save_flow_int)) then
            call save_flow(.false.)
        end if
    else
        if (mod(time_step, int(save_flow_int)) == 0) then
            call save_flow(.false.)
        end if
    end if

! Check for wall-time restriction
    call end_run_mpi(flag)
    if (flag) exit
    
end do

! Calculate and display the runtime for the simulation
call wall_time(end_time)
if (rank == 0) then
    write(*,*) 'Elapsed time (sec): ', end_time - start_time
    write(*,*) 'Seconds per iteration: ', (end_time - start_time)/n_time_steps
end if

time_step = time_step - 1
call save_flow(.true.)
if (save_stats_int < 0) then
    if ((mod(time, abs(save_stats_int)) + delta_t) <= abs(save_stats_int)) then
        call save_stats(.true.)
    end if
else
    if (mod(time_step,int(save_stats_int)) /= 0) then
        call save_stats(.true.)
    end if
end if

if (rank == 0) then
    write(*,*)
    write(*,*) '        ****** Hello world!  Have a nice day! ******'
    write(*,*)
end if
call mpi_barrier(mpi_comm_world,ierror)


end