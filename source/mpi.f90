

!----*|--.---------.---------.---------.---------.---------.---------.-|------
subroutine init_mpi
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header'

integer dims(2), perdim(2)
integer comm_cart
! This subroutine initializes all mpi variables

call mpi_init(ierror)

call mpi_comm_size(mpi_comm_world,nprocs,ierror)
call mpi_comm_rank(mpi_comm_world,rank,ierror)

np_s=nprocs**0.5
! rankz=mod(rank,np_s)
! ranky=(rank-rankz)/np_s

dims=(/np_s,np_s/)
perdim=(/1,1/)

call mpi_cart_create(mpi_comm_world, 2, dims, perdim, .false., &
                        comm_cart,ierror)
! In perdim, we put the information for the remain_dims
perdim=(/1,0/)
call mpi_cart_sub(comm_cart,perdim,mpi_comm_y,ierror)
perdim=(/0,1/)
call mpi_cart_sub(comm_cart,perdim,mpi_comm_z,ierror)

call mpi_comm_rank(mpi_comm_y,ranky,ierror)
call mpi_comm_rank(mpi_comm_z,rankz,ierror)

if (rank == 0) write(*,*) 'MPI initialized. nprocs: ',nprocs
call mpi_barrier(mpi_comm_world,ierror)
write(*,*) 'rank, ranky, rankz: ', rank, ranky, rankz

! Set a string to determine which input/output files to use
! When MPI is used, each process will read/write to files with the
! number of their rank (+1) appended to the end of the file.
! The string MPI_IO_NUM will be used to define the RANK+1 for each process
if (nprocs <= 10) then
    mpi_io_num=char(mod(rank+1,10)+48)
else if (nprocs <= 100) then
    mpi_io_num=char(mod(rank+1,100)/10+48) &
                //char(mod(rank+1,10)+48)
else if (nprocs <= 1000) then
    mpi_io_num=char(mod(rank+1,1000)/100+48) &
                //char(mod(rank+1,100)/10+48) &
                //char(mod(rank+1,10)+48)
else if (nprocs <= 10000) then
    mpi_io_num=char(mod(rank+1,10000)/1000+48) &
                //char(mod(rank+1,1000)/100+48) &
                //char(mod(rank+1,100)/10+48) &
                //char(mod(rank+1,10)+48)
else
    write(*,*) 'error, nprocs>10,000, unsupported problem size'
end if
mpi_num=rank

nxm_s=nx_s-1
nym_s=ny_s-1
nzm_s=nz_s-1

call mpi_barrier(mpi_comm_world,ierror)

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
subroutine transpose_mpi_z_to_x(ca_z,ca_x)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header'
! This subroutine initializes all mpi variables


complex*16 ca_z(0:nx_s/2,0:nz+1,0:ny_s)
complex*16 ca_x(0:nx/2,0:nz_s,0:ny_s)

complex*16 sendbuf(0:nkx_s,0:nz_s,0:ny_s)
complex*16 recvbuf(0:nkx_s,0:nz_s,0:ny_s)

integer diag, i, j, k
integer block
integer sendcount, recvcount
integer shift, dir, dir_min
integer i1, i2, j1, j2, k1, k2
integer ishift, jshift, kshift

! Zero Output array
ca_x=(0.d0,0.d0)
sendbuf=(0.d0,0.d0)
recvbuf=(0.d0,0.d0)

diag=mod(rank,np_s)

sendcount=(nkx_s+1)*(nz_s+1)*(ny_s+1)
recvcount=sendcount

do shift=0,np_s/2

    if ((shift*2 <= (np_s-1)).and.(shift /= 0)) then
        dir_min=-1
    else
        dir_min=1
    end if

    do dir=dir_min,1,2

        if (shift == 0) then
            block=diag
        else
            block=diag + shift*dir*(1-2*mod(diag/shift,2))
        end if

        if (block > np_s-1) then
            block=mod(block,np_s)
        else if (block < 0) then
            block=block+np_s
        end if

! Get a block from the source array to be sent to another process
        sendbuf=0.d0
        k1=(nz_s+1)*block
        k2=nz+1-k1
        do j=0,ny_s
            do k=0,min(k2,nz_s)
                do i=0,nkx_s
                    sendbuf(i,k,j)=ca_z(i,k+k1,j)
                end do
            end do
        end do

        if (block /= diag) then
! Communicate only if we aren't talking to ourself
            recvbuf=(0.d0,0.d0)
            call mpi_sendrecv( &
                sendbuf(0,0,0), sendcount, mpi_double_complex, &
                block+(rank/np_s)*np_s, 12, &
                recvbuf(0,0,0), recvcount, mpi_double_complex, &
                block+(rank/np_s)*np_s, 12, &
                mpi_comm_world, status, ierror)

! Place the newly recieved data into a temporary storage array
            i1 = (nkx_s+1)*block
            i2 = nkx_s+i1
            do j = 0,ny_s
                do k = 0,nz_s
                    do i = i1,min(i2,nx/2)
                        ca_x(i,k,j) = recvbuf(i-i1,k,j)
                    end do
                end do
            end do

        else
! Here block=diag, we just need to copy locally from A->B
            if (rank == 0) then
! If we are on RANK=0, just copy directly from input to output
                do j = 0,ny_s
                    do k = 0,nz_s
                        do i = 0,nkx_s
                            ca_x(i,k,j) = ca_z(i,k,j)
                        end do
                    end do
                end do
            else
! Here, RANK /= 0, so we need to shift the input and output
                k1 = (nz_s+1)*block
                k2 = nz+1-k1
                i1 = (nkx_s+1)*block
                i2 = i1+nx_s/2
                do j = 0,ny_s
                    do k = 0,min(k2,nz_s)
                        do i = i1,min(i2,nx/2)
                            ca_x(i,k,j) = ca_z(i-i1,k+k1,j)
                        end do
                    end do
                end do
            end if

        end if

    end do

end do

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
subroutine transpose_mpi_x_to_z(ca_x,ca_z)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header'
! This subroutine initializes all mpi variables

complex*16 ca_z(0:nkx_s,0:nz+1,0:ny_s)
complex*16 ca_x(0:nx/2,0:nz_s,0:ny_s)

complex*16 sendbuf(0:nkx_s,0:nz_s,0:ny_s)
complex*16 recvbuf(0:nkx_s,0:nz_s,0:ny_s)

integer diag,i,j,k
integer block
integer shift,dir,dir_min
integer sendcount,recvcount
integer i1,i2,j1,j2,k1,k2
integer ishift,jshift,kshift

! Zero Output array
ca_z=(0.d0,0.d0)
sendbuf=(0.d0,0.d0)
recvbuf=(0.d0,0.d0)

diag=mod(rank,np_s)

sendcount=(nkx_s+1)*(nz_s+1)*(ny_s+1)
recvcount=sendcount

do shift=0,np_s/2

    if ((shift*2 <= (np_s-1)).and.(shift /= 0)) then
        dir_min=-1
    else
        dir_min=1
    end if

    do dir=dir_min,1,2

        if (shift == 0) then
            block=diag
        else
            block=diag + shift*dir*(1-2*mod(diag/shift,2))
        end if

        if (block > np_s-1) then
            block=mod(block,np_s)
        else if (block < 0) then
            block=block+np_s
        end if

! Get a block from the source array to be sent to another process
        sendbuf = 0.d0
        i1 = (nkx_s+1)*block
        i2 = nx/2-i1
        do j = 0,ny_s
            do k = 0,nz_s
                do i = 0,min(i2,nkx_s)
                    sendbuf(i,k,j) = ca_x(i+i1,k,j)
                end do
            end do
        end do

        if (block /= diag) then
! Communicate only if we aren't talking to ourself
            recvbuf = (0.d0,0.d0)
            call mpi_sendrecv( &
                sendbuf(0,0,0), sendcount, mpi_double_complex, &
                block+(rank/np_s)*np_s, 13, &
                recvbuf(0,0,0), recvcount, mpi_double_complex, &
                block+(rank/np_s)*np_s, 13, &
                mpi_comm_world, status, ierror)

! Place the newly recieved data into a temporary storage array
            k1 = (nz_s+1)*block
            k2 = nz_s+k1
            do j = 0,ny_s
                do k = k1,min(k2,nz+1)
                    do i = 0,nkx_s
! Make sure that we don't write beyond the end of the array
                        ca_z(i,k,j) = recvbuf(i,k-k1,j)
                    end do
                end do
            end do

        else
! Here block = diag, we just need to copy locally from A->B
            if (rank == 0) then
                do j = 0,ny_s
                    do k = 0,nz_s
                        do i = 0,nkx_s
                            ca_z(i,k,j) = ca_x(i,k,j)
                        end do
                    end do
                end do
            else
                i1 = (nkx_s+1)*block
                k1 = (nz_s+1)*block
                k2 = nz_s+k1
                do j = 0,ny_s
                    do k = k1,min(k2,nz+1)
                        do i = 0,nkx_s
                            ca_z(i,k,j) = ca_x(i+i1,k-k1,j)
                        end do
                    end do
                end do
            end if

        end if

    end do

end do

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
subroutine transpose_mpi_x_to_y(a_x,a_y,nx,ny,nz,nx_s,ny_s,nz_s)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header_mpi'
! This subroutine initializes all mpi variables

integer nx,nx_s,ny,ny_s,nz,nz_s
complex*16 a_x(0:nx,0:nz_s,0:ny_s)
complex*16 a_y(0:nx_s,0:nz_s,0:ny)
complex*16 b(0:nx_s,0:nz_s,0:ny)

complex*16 sendbuf(0:nx_s,0:nz_s,0:ny_s)
complex*16 recvbuf(0:nx_s,0:nz_s,0:ny_s)

integer diag,i,j,k
integer block
integer sendcount,recvcount

diag = mod(rank,np_s)

sendcount = (nx_s+1)*(nz_s+1)*(ny_s+1)
recvcount = sendcount

do block = 0,np_s-1

! Get a block from the source array to be sent to another process
    sendbuf = 0.d0
    do j = 0,ny_s
        do k = 0,nz_s
            do i = 0,nx_s
                sendbuf(i,k,j) = a_x(i+(nx_s+1)*block,k,j)
            end do
        end do
    end do

    if (block /= diag) then
! Communicate only if we aren't talking to ourself
        recvbuf = (0.d0,0.d0)
        call mpi_sendrecv( &
            sendbuf,sendcount,mpi_double_complex, &
            rank+block-diag,rank+block, &
            recvbuf,recvcount,mpi_double_complex, &
            rank+block-diag,rank+block, &
            mpi_comm_world,status,ierror)

! Place the newly recieved data into a temporary storage array
        do i = 0,nx_s
            do k = 0,nz_s
                do j = 0,ny_s
                    if (j+(ny_s+1)*block <= ny) then
! Make sure that we don't write beyond the end of the array
                        b(i,k,j+(ny_s+1)*block) = recvbuf(i,k,j)
                    end if
                end do
            end do
        end do

    else
! Here block=diag, we just need to copy locally from A->B

        do i = 0,nx_s
            do k = 0,nz_s
                do j = 0,ny_s
                    sendbuf(i,k,j) = a_x(i+(nx_s+1)*block,k,j)
                end do
            end do
        end do

        do i = 0,nx_s
            do k = 0,nz_s
                do j = 0,ny_s
                    if (j+(ny_s+1)*block <= ny) then
                        b(i,k,j+(ny_s+1)*block) = sendbuf(i,k,j)
                    end if
                end do
            end do
        end do
    end if

end do

! Finally, copy the temporary storage array to the final array
do i = 0,nx_s
    do k = 0,nz_s
        do j = 0,ny
            a_y(i,k,j) = b(i,k,j)
        end do
    end do
end do

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
subroutine transpose_mpi_y_to_x(a_y,a_x,nx,ny,nz,nx_s,ny_s,nz_s)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header_mpi'
! This subroutine initializes all mpi variables

integer nx,nx_s,ny,ny_s,nz,nz_s
complex*16 a_x(0:nx,0:nz_s,0:ny_s)
complex*16 a_y(0:nx_s,0:nz_s,0:ny)
complex*16 b(0:nx,0:nz_s,0:ny_s)

complex*16 sendbuf(0:nx_s,0:nz_s,0:ny_s)
complex*16 recvbuf(0:nx_s,0:nz_s,0:ny_s)

integer diag,i,j,k
integer block
integer sendcount,recvcount

integer i1,i2,j1,j2,k1,k2

diag = mod(rank,np_s)

sendcount = (nx_s+1)*(nz_s+1)*(ny_s+1)
recvcount = sendcount

do block = 0,np_s-1

! Get a block from the source array to be sent to another process
    sendbuf = 0.d0
    do i = 0,nx_s
        do k = 0,nz_s
            do j = 0,ny_s
                sendbuf(i,k,j) = a_y(i,k,j+(ny_s+1)*block)
            end do
        end do
    end do

    if (block /= diag) then
        recvbuf = (0.d0,0.d0)
! Communicate only if we aren't talking to ourself
        call mpi_sendrecv( &
            sendbuf, sendcount, mpi_double_complex, &
            rank+block-diag, rank+block, &
            recvbuf, recvcount, mpi_double_complex, &
            rank+block-diag, rank+block, &
            mpi_comm_world, status, ierror)

! Place the newly recieved data into a temporary storage array
        i2 = nx-block*(nx_s+1)
        do i = 0,min(i2,nx_s)
            do k = 0,nz_s
                do j = 0,ny_s
! Make sure that we don't write beyond the end of the array
                    b(i+(nx_s+1)*block,k,j) = recvbuf(i,k,j)
                end do
            end do
        end do

    else
! Here block=diag, we just need to copy locally from A->B

        do i = 0,nx_s
            do k = 0,nz_s
                do j = 0,ny_s
                    sendbuf(i,k,j) = a_y(i,k,j+(ny_s+1)*block)
                end do
            end do
        end do
        do i = 0,nx_s
            do k = 0,nz_s
                do j = 0,ny_s
                    if (i+(nx_s+1)*block <= nx) then
                        b(i+(nx_s+1)*block,k,j) = sendbuf(i,k,j)
                    end if
                end do
            end do
        end do
    end if

end do

! Finally, copy the temporary storage array to the final array
do i = 0,nx
    do k = 0,nz_s
        do j = 0,ny_s
            a_x(i,k,j) = b(i,k,j)
        end do
    end do
end do

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
subroutine transpose_mpi_y_to_z(ca_y,ca_z)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header'
! This subroutine initializes all mpi variables


complex*16 ca_z(0:nx_s/2,0:nz+1,0:ny_s)
complex*16 ca_y(0:nx_s/2,0:nz_s,0:ny+1)

complex*16 sendbuf(0:nkx_s,0:tnkz_s,0:ny_s)
complex*16 recvbuf(0:nkx_s,0:tnkz_s,0:ny_s)

integer diag,i,j,k
integer block
integer sendcount,recvcount
integer shift,dir,dir_min
integer i1,i2,k1,k2,j1,j2
integer ishift,jshift,kshift

! Zero Output array
ca_z = (0.d0,0.d0)
sendbuf = (0.d0,0.d0)
recvbuf = (0.d0,0.d0)

diag = rank/np_s

sendcount = (nkx_s+1)*(tnkz_s+1)*(ny_s+1)
recvcount = sendcount

do shift = 0,np_s/2


    if ((shift*2 <= (np_s-1)).and.(shift /= 0)) then
        dir_min = -1
    else
        dir_min = 1
    end if

    do dir = dir_min,1,2

        if (shift == 0) then
            block = diag
        else
            block = diag + shift*dir*(1-2*mod(diag/shift,2))
        end if

        if (block > np_s-1) then
            block = mod(block,np_s)
        else if (block < 0) then
            block = block+np_s
        end if

! Get a block from the source array to be sent to another process
        sendbuf = 0.d0
        j1 = (ny_s+1)*block
        j2 = ny+1-j1
        do j = 0,min(j2,ny_s)
            do k = 0,tnkz_s
                do i = 0,nkx_s
                    sendbuf(i,k,j) = ca_y(i,k,j+j1)
                end do
            end do
        end do


        if (block /= diag) then
! Communicate only if we aren't talking to ourself
            recvbuf = (0.d0,0.d0)
            call mpi_sendrecv( &
                    sendbuf,sendcount,mpi_double_complex, &
                    mod(rank,np_s)+block*np_s, &
                    rank+mod(rank,np_s)+block*np_s, &
                    recvbuf,recvcount,mpi_double_complex, &
                    mod(rank,np_s)+block*np_s, &
                    rank+mod(rank,np_s)+block*np_s, &
                    mpi_comm_world,status,ierror)

! Place the newly recieved data into a temporary storage array
            k1 = (tnkz_s+1)*block
            k2 = tnkz_s+k1
            do j = 0,ny_s
                do k = k1,min(k2,nz+1)
                    do i = 0,nkx_s
                        ca_z(i,k,j) = recvbuf(i,k-k1,j)
                    end do
                end do
            end do

        else
! Here block=diag, we just need to copy locally from A->B
            if (rank == 0) then
                do j = 0,ny_s
                    do k = 0,tnkz_s
                        do i = 0,nkx_s
                            ca_z(i,k,j) = ca_y(i,k,j)
                        end do
                    end do
                end do
            else
                j1 = (ny_s+1)*block
                j2 = ny+1-j1
                k1 = (tnkz_s+1)*block
                k2 = k1+tnkz_s
                do j = 0,min(j2,ny_s)
                    do k = k1,min(k2,nz+1)
                        do i = 0,nkx_s
                            ca_z(i,k,j) = ca_y(i,k-k1,j+j1)
                        end do
                    end do
                end do
            end if

        end if


    end do

end do

return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
subroutine transpose_mpi_z_to_y(ca_z,ca_y)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header'
! This subroutine initializes all mpi variables

complex*16 ca_z(0:nkx_s,0:nz+1,0:ny_s)
complex*16 ca_y(0:nkx_s,0:tnkz_s,0:ny+1)

complex*16 sendbuf(0:nkx_s,0:tnkz_s,0:ny_s)
complex*16 recvbuf(0:nkx_s,0:tnkz_s,0:ny_s)

integer diag,i,j,k
integer block
integer sendcount,recvcount
integer shift,dir,dir_min
integer i1,i2,j1,j2,k1,k2
integer ishift,jshift,kshift

! Zero Output array, Working arrays
ca_y = (0.d0,0.d0)
sendbuf = (0.d0,0.d0)
recvbuf = (0.d0,0.d0)

diag = rank/np_s

sendcount = (nkx_s+1)*(tnkz_s+1)*(ny_s+1)
recvcount = sendcount

do shift = 0,np_s/2

    if ((shift*2 <= (np_s-1)).and.(shift /= 0)) then
        dir_min = -1
    else
        dir_min = 1
    end if

    do dir = dir_min,1,2

        if (shift == 0) then
            block = diag
        else
            block = diag + shift*dir*(1-2*mod(diag/shift,2))
        end if

        if (block > np_s-1) then
            block = mod(block,np_s)
        else if (block < 0) then
            block = block+np_s
        end if

! Get a block from the source array to be sent to another process
        sendbuf = 0.d0
        k1 = (tnkz_s+1)*block
        k2 = nz+1-k1
        do j = 0,ny_s
            do k = 0,min(k2,tnkz_s)
                do i = 0,nkx_s
                    sendbuf(i,k,j) = ca_z(i,k+k1,j)
                end do
            end do
        end do

        if (block /= diag) then
            recvbuf = (0.d0,0.d0)
! Communicate only if we aren't talking to ourself
            call mpi_sendrecv( &
                sendbuf, sendcount, mpi_double_complex, &
                mod(rank,np_s)+block*np_s, 14, &
                recvbuf, recvcount, mpi_double_complex, &
                mod(rank,np_s)+block*np_s, 14, &
                mpi_comm_world, status, ierror)

! Place the newly recieved data into a temporary storage array
            j1 = (ny_s+1)*block
            j2 = ny_s+j1
            do j = j1,min(j2,ny+1)
                do k = 0,tnkz_s
                    do i = 0,nkx_s
! Make sure that we don't write beyond the end of the array
                        ca_y(i,k,j) = recvbuf(i,k,j-j1)
                    end do
                end do
            end do

        else
! Here block=diag, we just need to copy locally from A->B
            if (rank == 0) then
                do j = 0,ny_s
                    do k = 0,tnkz_s
                        do i = 0,nkx_s
                            ca_y(i,k,j) = ca_z(i,k,j)
                        end do
                    end do
                end do
            else
                k1 = (tnkz_s+1)*block
                k2 = nz+1-k1
                j1 = (ny_s+1)*block
                j2 = ny_s+j1
                do j = j1,min(j2,ny+1)
                    do k = 0,min(k2,tnkz_s)
                        do i = 0,nkx_s
                            ca_y(i,k,j) = ca_z(i,k+k1,j-j1)
                        end do
                    end do
                end do
            end if
        end if

    end do

end do


return
end

!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
subroutine courant_mpi
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
! This subroutine sets the timestep based on the specified CFL number
! The subroutine should be called with the velocity in physical space

include 'header'

real*8 vel
real*8 dt_local,dt_test
real*8 dt_x,dt_y,dt_z
real*8 dt
integer i,j,k
integer imin,jmin,kmin

! Set the initial dt to some arbitrary large number
! For chemotaxis with flow
dt_local = 7.d-3
! For no flow: (rescaled for 128^3)
! dt_local = 7.d-2*25.d0
! Set the maximum dt based on the buoyancy period
! dt_local = 0.005d0*2.d0*PI/sqrt(RI_TAU(1))

dt_test = 0.d0

do j = 0,ny_s
    do k = 0,nz_s
        do i = 0,nxm
            if (u1(i,k,j) /= 0.d0) then
                dt_x = abs(cfl*dx_th(i)/abs(u1(i,k,j)))
            else
                dt_x = 999.d0
            end if
            if (u2(i,k,j) /= 0.d0) then
                dt_y = abs(cfl*dy_th(j)/abs(u2(i,k,j)))
            else
                dt_y = 999.d0
            end if
! Note, thermal wind advection included
!           dt_z = abs(cfl*dz(k)/abs( &
!                   U3(i,k,j)+(RI_TAU(N)/I_RO_TAU)*DRHODX*GYF(J) &
!                   ))
            if (u3(i,k,j) /= 0.d0) then
                dt_z = abs(cfl*dz_th(k)/abs(u3(i,k,j)))
            else
                dt_z = 999.d0
            end if
            dt_local = min(dt_local,dt_x,dt_y,dt_z)
        end do
    end do
end do

! Now we have the minimum DELTA_T on each process
! share the information to obtain the global minimum

call mpi_allreduce(dt_local, dt, 1, &
        mpi_double_precision, mpi_min, mpi_comm_world, ierror)

if (dt <= 0) then
    write(*,*) 'Error: dt<=0 in courant'
! Set DELTA_T to some small default value
    delta_t = 0.0001d0
else if (dt >= 1000.) then
    write(*,*) 'warning: delta_t > 1000, value capped at 1000'
    delta_t = 1000.d0
    h_bar(1) = delta_t*(8.0/15.0)
    h_bar(2) = delta_t*(2.0/15.0)
    h_bar(3) = delta_t*(5.0/15.0)
else
    delta_t = dt
    h_bar(1) = delta_t*(8.0/15.0)
    h_bar(2) = delta_t*(2.0/15.0)
    h_bar(3) = delta_t*(5.0/15.0)
end if

return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
subroutine end_run_mpi(flag)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header'

logical flag

if (rank == 0) then
    call end_run(flag)
end if
call mpi_bcast(flag, 1, mpi_logical, 0, mpi_comm_world, ierror)

end