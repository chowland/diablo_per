!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
subroutine init_mpi_th
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header'
! This subroutine initializes all mpi variables

nxm_s_th = nx_s_th-1
nym_s_th = ny_s_th-1
nzm_s_th = nz_s_th-1

return
end


!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
subroutine transpose_mpi_z_to_x_th(ca_z,ca_x)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header'
! This subroutine initializes all mpi variables


complex*16 ca_z(0:nx_s_th/2,0:nz_th+1,0:ny_s_th)
complex*16 ca_x(0:nx_th/2,0:nz_s_th,0:ny_s_th)

complex*16 sendbuf(0:nkx_s_th,0:nz_s_th,0:ny_s_th)
complex*16 recvbuf(0:nkx_s_th,0:nz_s_th,0:ny_s_th)

integer diag,i,j,k
integer block
integer sendcount,recvcount
integer shift,dir,dir_min
integer i1,i2,j1,j2,k1,k2
integer ishift,jshift,kshift

! Zero Output array
ca_x = (0.d0,0.d0)
sendbuf = (0.d0,0.d0)
recvbuf = (0.d0,0.d0)


diag = mod(rank,np_s)

sendcount = (nkx_s_th+1)*(nz_s_th+1)*(ny_s_th+1)
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
        k1 = (nz_s_th+1)*block
        k2 = nz_th+1-k1
        do j = 0,ny_s_th
            do k = 0,min(k2,nz_s_th)
                do i = 0,nkx_s_th
                    sendbuf(i,k,j) = ca_z(i,k+k1,j)
                end do
            end do
        end do

        if (block /= diag) then
! Communicate only if we aren't talking to ourself
            recvbuf = 0.d0
            call mpi_sendrecv( &
                    sendbuf(0,0,0),sendcount,mpi_double_complex, &
                    block+(rank/np_s)*np_s,12, &
                    recvbuf(0,0,0),recvcount,mpi_double_complex, &
                    block+(rank/np_s)*np_s,12, &
                    mpi_comm_world,status,ierror)

! Place the newly recieved data into a temporary storage array
            i1 = (nkx_s_th+1)*block
            i2 = nkx_s_th+i1
            do j = 0,ny_s_th
                do k = 0,nz_s_th
                    do i = i1,min(i2,nx_th/2)
! Make sure that we don't write beyond the end of the array
                        ca_x(i,k,j) = recvbuf(i-i1,k,j)
                    end do
                end do
            end do

        else
! Here block = diag, we just need to copy locally from A->B
            if (rank == 0) then
! If we are on RANK = 0, just copy directly from input to output
                do j = 0,ny_s_th
                    do k = 0,nz_s_th
                        do i = 0,nkx_s_th
                            ca_x(i,k,j) = ca_z(i,k,j)
                        end do
                    end do
                end do
            else
! Here, RANK /= 0, so we need to shift the input and output
                k1 = (nz_s_th+1)*block
                k2 = nz_th+1-k1
                i1 = (nkx_s_th+1)*block
                i2 = i1+nkx_s_th
                do j = 0,ny_s_th
                    do k = 0,min(k2,nz_s_th)
                        do i = i1,min(i2,nx_th/2)
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
subroutine transpose_mpi_x_to_z_th(ca_x,ca_z)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header'
! This subroutine initializes all mpi variables

complex*16 ca_z(0:nkx_s_th,0:nz_th+1,0:ny_s_th)
complex*16 ca_x(0:nx_th/2,0:nz_s_th,0:ny_s_th)

complex*16 sendbuf(0:nkx_s_th,0:nz_s_th,0:ny_s_th)
complex*16 recvbuf(0:nkx_s_th,0:nz_s_th,0:ny_s_th)

integer diag,i,j,k
integer block
integer shift,dir,dir_min
integer sendcount,recvcount
integer i1,i2,j1,j2,k1,k2
integer ishift,jshift,kshift

! Zero Output array
ca_z = (0.d0,0.d0)
sendbuf = (0.d0,0.d0)
recvbuf = (0.d0,0.d0)

diag = mod(rank,np_s)

sendcount = (nkx_s_th+1)*(nz_s_th+1)*(ny_s_th+1)
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
        i1 = (nkx_s_th+1)*block
        i2 = nx_th/2-i1
        do j = 0,ny_s_th
            do k = 0,nz_s_th
                do i = 0,min(i2,nkx_s_th)
                    sendbuf(i,k,j) = ca_x(i+i1,k,j)
                end do
            end do
        end do

        if (block /= diag) then
! Communicate only if we aren't talking to ourself
            recvbuf = 0.d0
            call mpi_sendrecv( &
                sendbuf(0,0,0),sendcount,mpi_double_complex, &
                block+(rank/np_s)*np_s,13, &
                recvbuf(0,0,0),recvcount,mpi_double_complex, &
                block+(rank/np_s)*np_s,13, &
                mpi_comm_world,status,ierror)

! Place the newly recieved data into a temporary storage array
            k1 = (nz_s_th+1)*block
            k2 = nz_s_th+k1
            do j = 0,ny_s_th
                do k = k1,min(k2,nz_th+1)
                    do i = 0,nkx_s_th
! Make sure that we don't write beyond the end of the array
                        ca_z(i,k,j) = recvbuf(i,k-k1,j)
                    end do
                end do
            end do

        else
! Here block = diag, we just need to copy locally from A->B
            if (rank == 0) then
                do j = 0,ny_s_th
                    do k = 0,nz_s_th
                        do i = 0,nkx_s_th
                            ca_z(i,k,j) = ca_x(i,k,j)
                        end do
                    end do
                end do
            else
                i1 = (nkx_s_th+1)*block
                k1 = (nz_s_th+1)*block
                k2 = nz_s_th+k1
                do j = 0,ny_s_th
                    do k = k1,min(k2,nz_th+1)
                        do i = 0,nkx_s_th
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
subroutine transpose_mpi_y_to_z_th(ca_y,ca_z)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header'
! This subroutine initializes all mpi variables


complex*16 ca_z(0:nx_s_th/2,0:nz_th+1,0:ny_s_th)
complex*16 ca_y(0:nx_s_th/2,0:nz_s_th,0:ny_th+1)

complex*16 sendbuf(0:nkx_s_th,0:tnkz_s_th,0:ny_s_th)
complex*16 recvbuf(0:nkx_s_th,0:tnkz_s_th,0:ny_s_th)

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

sendcount = (nkx_s_th+1)*(tnkz_s_th+1)*(ny_s_th+1)
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
        j1 = (ny_s_th+1)*block
        j2 = ny_th+1-j1
        do j = 0,min(j2,ny_s_th)
            do k = 0,tnkz_s_th
                do i = 0,nkx_s_th
                    sendbuf(i,k,j) = ca_y(i,k,j+j1)
                end do
            end do
        end do

        if (block /= diag) then
! Communicate only if we aren't talking to ourself
            recvbuf = (0.d0,0.d0)
            call mpi_sendrecv( &
                sendbuf,sendcount,mpi_double_complex, &
                mod(rank,np_s)+block*np_s,15, &
                recvbuf,recvcount,mpi_double_complex, &
                mod(rank,np_s)+block*np_s,15, &
                mpi_comm_world,status,ierror)

! Place the newly recieved data into a temporary storage array
            k1 = (tnkz_s_th+1)*block
            k2 = tnkz_s_th+k1
            do j = 0,ny_s_th
                do k = k1,min(k2,nz_th+1)
                    do i = 0,nkx_s_th
! Make sure that we don't write beyond the end of the array
                        ca_z(i,k,j) = recvbuf(i,k-k1,j)
                    end do
                end do
            end do

        else
! Here block = diag, we just need to copy locally from A->B
            if (rank == 0) then
                do j = 0,ny_s_th
                    do k = 0,tnkz_s_th
                        do i = 0,nkx_s_th
                            ca_z(i,k,j) = ca_y(i,k,j)
                        end do
                    end do
                end do
            else
                j1 = (ny_s_th+1)*block
                j2 = ny_th+1-j1
                k1 = (tnkz_s_th+1)*block
                k2 = k1+tnkz_s_th
                do j = 0,min(j2,ny_s_th)
                    do k = k1,min(k2,nz_th+1)
                        do i = 0,nkx_s_th
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
subroutine transpose_mpi_z_to_y_th(ca_z,ca_y)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header'
! This subroutine initializes all mpi variables

complex*16 ca_z(0:nkx_s_th,0:nz_th+1,0:ny_s_th)
complex*16 ca_y(0:nkx_s_th,0:tnkz_s_th,0:ny_th+1)

complex*16 sendbuf(0:nkx_s_th,0:tnkz_s_th,0:ny_s_th)
complex*16 recvbuf(0:nkx_s_th,0:tnkz_s_th,0:ny_s_th)

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

sendcount = (nkx_s_th+1)*(tnkz_s_th+1)*(ny_s_th+1)
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
        k1 = (tnkz_s_th+1)*block
        k2 = nz_th+1-k1
        do j = 0,ny_s_th
            do k = 0,min(k2,tnkz_s_th)
                do i = 0,nkx_s_th
                    sendbuf(i,k,j) = ca_z(i,k+k1,j)
                end do
            end do
        end do

        if (block /= diag) then
! Communicate only if we aren't talking to ourself
            recvbuf = 0.d0
            call mpi_sendrecv( &
                sendbuf,sendcount,mpi_double_complex, &
                mod(rank,np_s)+block*np_s,14, &
                recvbuf,recvcount,mpi_double_complex, &
                mod(rank,np_s)+block*np_s,14, &
                mpi_comm_world,status,ierror)

! Place the newly recieved data into a temporary storage array
            j1 = (ny_s_th+1)*block
            j2 = ny_s_th+j1
            do j = j1,min(j2,ny_th+1)
                do k = 0,tnkz_s_th
                    do i = 0,nkx_s_th
! Make sure that we don't write beyond the end of the array
                        ca_y(i,k,j) = recvbuf(i,k,j-j1)
                    end do
                end do
            end do

        else
! Here block = diag, we just need to copy locally from A->B
            if (rank == 0) then
                do j = 0,ny_s_th
                    do k = 0,tnkz_s_th
                        do i = 0,nkx_s_th
                            ca_y(i,k,j) = ca_z(i,k,j)
                        end do
                    end do
                end do
            else
                k1 = (tnkz_s_th+1)*block
                k2 = nz_th+1-k1
                j1 = (ny_s_th+1)*block
                j2 = ny_s_th+j1
                do j = j1,min(j2,ny_th+1)
                    do k = 0,min(k2,tnkz_s_th)
                        do i = 0,nkx_s_th
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



! The following subroutine is for taking a velocity field in Fourier space
! and interpolating it to the scalar (TH) grid in Physical space
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
subroutine transpose_mpi_y_to_z_interp(ca_y,ca_z)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header'
! This subroutine initializes all mpi variables


complex*16 ca_z(0:nx_s/2,0:nz+1,0:ny_s_th)
complex*16 ca_y(0:nx_s/2,0:nz_s,0:ny_th+1)

complex*16 sendbuf(0:nkx_s,0:tnkz_s,0:ny_s_th)
complex*16 recvbuf(0:nkx_s,0:tnkz_s,0:ny_s_th)

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

sendcount = (nkx_s+1)*(tnkz_s+1)*(ny_s_th+1)
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
        j1 = (ny_s_th+1)*block
        j2 = ny_th+1-j1
        do j = 0,min(j2,ny_s_th)
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
                mod(rank,np_s)+block*np_s,15, &
                recvbuf,recvcount,mpi_double_complex, &
                mod(rank,np_s)+block*np_s,15, &
                mpi_comm_world,status,ierror)

! Place the newly recieved data into a temporary storage array
            k1 = (tnkz_s+1)*block
            k2 = tnkz_s+k1
            do j = 0,ny_s_th
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
                do j = 0,ny_s_th
                    do k = 0,tnkz_s
                        do i = 0,nkx_s
                            ca_z(i,k,j) = ca_y(i,k,j)
                        end do
                    end do
                end do
            else
                j1 = (ny_s_th+1)*block
                j2 = ny_th+1-j1
                k1 = (tnkz_s+1)*block
                k2 = k1+tnkz_s
                do j = 0,min(j2,ny_s_th)
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
subroutine transpose_mpi_z_to_x_interp(ca_z,ca_x)
!----*|--.---------.---------.---------.---------.---------.---------.-|-----|
include 'header'
! This subroutine initializes all mpi variables


complex*16 ca_z(0:nx_s/2,0:nz_th+1,0:ny_s_th)
complex*16 ca_x(0:nx/2,0:nz_s_th,0:ny_s_th)

complex*16 sendbuf(0:nkx_s,0:nz_s_th,0:ny_s_th)
complex*16 recvbuf(0:nkx_s,0:nz_s_th,0:ny_s_th)

integer diag,i,j,k
integer block
integer sendcount,recvcount
integer shift,dir,dir_min
integer i1,i2,j1,j2,k1,k2
integer ishift,jshift,kshift

! Zero Output array
ca_x = (0.d0,0.d0)
sendbuf = (0.d0,0.d0)
recvbuf = (0.d0,0.d0)


diag = mod(rank,np_s)

sendcount = (nkx_s+1)*(nz_s_th+1)*(ny_s_th+1)
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
        k1 = (nz_s_th+1)*block
        k2 = nz_th+1-k1
        do j = 0,ny_s_th
            do k = 0,min(k2,nz_s_th)
                do i = 0,nkx_s
                    sendbuf(i,k,j) = ca_z(i,k+k1,j)
                end do
            end do
        end do

        if (block /= diag) then
! Communicate only if we aren't talking to ourself
            recvbuf = 0.d0
            call mpi_sendrecv( &
                sendbuf(0,0,0),sendcount,mpi_double_complex, &
                block+(rank/np_s)*np_s,12, &
                recvbuf(0,0,0),recvcount,mpi_double_complex, &
                block+(rank/np_s)*np_s,12, &
                mpi_comm_world,status,ierror)

! Place the newly recieved data into a temporary storage array
            i1 = (nkx_s+1)*block
            i2 = nkx_s+i1
            do j = 0,ny_s_th
                do k = 0,nz_s_th
                    do i = i1,min(i2,nx/2)
! Make sure that we don't write beyond the end of the array
                        ca_x(i,k,j) = recvbuf(i-i1,k,j)
                    end do
                end do
            end do

        else
! Here block = diag, we just need to copy locally from A->B
            if (rank == 0) then
! If we are on RANK = 0, just copy directly from input to output
                do j = 0,ny_s_th
                    do k = 0,nz_s_th
                        do i = 0,nkx_s
                            ca_x(i,k,j) = ca_z(i,k,j)
                        end do
                    end do
                end do
            else
! Here, RANK /= 0, so we need to shift the input and output
                k1 = (nz_s_th+1)*block
                k2 = nz_th+1-k1
                i1 = (nkx_s+1)*block
                i2 = i1+nkx_s
                do j = 0,ny_s_th
                    do k = 0,min(k2,nz_s_th)
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