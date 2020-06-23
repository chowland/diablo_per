subroutine user_rhs_per_fourier

include 'header'
! Here, you can add terms to the right hand side of the momentum
! and scalar equations.
! The right hand side foring arrays, CF1, CF2, CF3, CFTH are in Fourier space.
! CF1, CF2, CF3 are available as working variables, but also define the output.
! CFTH already contains the background stratification, so is not a working variable.
! The velocity and scalars are available in Fourier space.
! CS1 and CSTH1 are available as working variables.

integer i,j,k,n

real*8 alpha, k0, f0, b2, puf, pff, puf_sum, pf, pf_sum
real*8 P_aim, f01, f02, kappa2, k2, Ps, ps_sum, fs1, fs2

P_aim=ri_tau(1)*target_reb*nu
Puf=0.d0
Ps=0.d0
Pf=0.d0
if (first_time) pff=0.d0
call random_seed

if (f_type.eq.1) then
! **** LOW WAVENUMBER (FURUE, 2003)-like FORCING ****

    k0=7.

    do j=0,tnky
        do k=0,tnkz_s
            do i=0,nkx_s
                if ((kx_s(i).eq.0) .and. ((kz_s(k).lt.0) .or. &
                    ((kz_s(k).eq.0) .and. (ky(j).lt.0)))) then
                else
                    if ((i+rankz*(nkx_s+1).le.nkx) .and. &
                            (k+ranky*(tnkz_s+1).le.tnkz)) then
                        kappa2=kx2_s(i)+kz2_s(k)
                        k2=kappa2+ky2(j)
    !                   if ((k2.le.100.) .and. (ky(j).ne.0)
                        if ((k2.gt.6.25) .and. (k2.lt. 12.25) &
                                .and. (kappa2.ne.0)) then
                            call random_number(alpha)
                            alpha=2.d0*pi*alpha ! random phase of forcing
                            cs1(i,k,j)=cexp(cmplx(0,alpha))
    !     &                        *k2**(1./4.)*(k2+k0**2)**(-3)
                            cf1(i,k,j)=cs1(i,k,j)*ky(j)*kx_s(i)/sqrt(k2*kappa2)
                            cf2(i,k,j)=-cs1(i,k,j)*sqrt(kappa2)/sqrt(k2)
                            cf3(i,k,j)=cs1(i,k,j)*ky(j)*kz_s(k)/sqrt(k2*kappa2)
                            csth1(i,k,j)=cs1(i,k,j)*ci/sqrt(ri_tau(1))
                            puf=puf+conjg(cu1(i,k,j))*cf1(i,k,j)+conjg(cu2(i,k,j)) &
                                    *cf2(i,k,j)+conjg(cu3(i,k,j))*cf3(i,k,j) &
                                    +ri_tau(1)*conjg(cth(i,k,j,1))*csth1(i,k,j)
                            if (first_time) then
                                pff=pff+1.d0 !abs(cs1(i,k,j))**2
                            end if
                        end if
                    end if
                end if
            end do
        end do
    end do

    call mpi_allreduce(puf,puf_sum,1,mpi_double_precision,mpi_sum, &
                        mpi_comm_world,ierror)
    if (first_time) then
        call mpi_allreduce(pff,pff_sum,1,mpi_double_precision,mpi_sum, &
                            mpi_comm_world,ierror)
    end if
    puf_sum=2.d0*puf_sum
    if (first_time) pff_sum=2.d0*pff_sum

    f01=(-puf_sum+sqrt(puf_sum**2+4*pff_sum*delta_t*p_aim)) &
            /(2*delta_t*pff_sum)
    f02=(-puf_sum-sqrt(puf_sum**2+4*pff_sum*delta_t*p_aim)) &
            /(2*delta_t*pff_sum)
    if (abs(f01)<abs(f02)) then
        do j=0,tnky
            do k=0,tnkz_s
                do i=0,nkx_s
                    cf1(i,k,j)=f01*cf1(i,k,j)
                    cf2(i,k,j)=f01*cf2(i,k,j)
                    cf3(i,k,j)=f01*cf3(i,k,j)
                    cfth(i,k,j,1)=cfth(i,k,j,1)+f01*csth1(i,k,j)
                end do
            end do
        end do
    else
        do j=0,tnky
            do k=0,tnkz_s
                do i=0,nkx_s
                    cf1(i,k,j)=f02*cf1(i,k,j)
                    cf2(i,k,j)=f02*cf2(i,k,j)
                    cf3(i,k,j)=f02*cf3(i,k,j)
                    cfth(i,k,j,1)=cfth(i,k,j,1)+f02*csth1(i,k,j)
                end do
            end do
        end do
    end if


else if (f_type.eq.2) then
! **** MAFFIOLI (2017) FORCING WITH CONSTANT POWER INJECTION ****

    do j=0,tnky
        do k=0,tnkz_s
            do i=0,nkx_s
                if ((kx_s(i).eq.0) .and. ((kz_s(k).lt.0) .or. &
                        ((kz_s(k).eq.0) .and. (ky(j).lt.0)))) then
                else
                    if ((i+rankz*(nkx_s+1).le.nkx) .and. &
                            (k+ranky*(tnkz_s+1).le.tnkz)) then
                        kappa2=kx2_s(i)+kz2_s(k)
                        if ((kappa2.gt.6.25) .and. (ky(j).eq.0) &
                                .and. (kappa2.lt.12.25)) then
                            call random_number(alpha)
                            cf1(i,k,j)=kz_s(k)*cexp(cmplx(0,2.d0*pi*alpha)) &
                                        /sqrt(pi)/3**1.5
                            cf3(i,k,j)=-kx_s(i)*cexp(cmplx(0,2.d0*pi*alpha)) &
                                        /sqrt(pi)/3**1.5
                            puf=puf+conjg(cu1(i,k,j))*cf1(i,k,j) &
                                    +conjg(cu3(i,k,j))*cf3(i,k,j)
                            if (first_time) then
                                pff=pff+0.5*(abs(cf1(i,k,j))**2+abs(cf3(i,k,j))**2)
                            end if
                        end if
                    end if
                end if
            end do
        end do
    end do

    call mpi_allreduce(puf,puf_sum,1,mpi_double_precision,mpi_sum, &
                        mpi_comm_world,ierror)
    if (first_time) then
        call mpi_allreduce(pff,pff_sum,1,mpi_double_precision,mpi_sum, &
                            mpi_comm_world,ierror)
    end if
    puf_sum=2.d0*puf_sum
    if (first_time) pff_sum=2.d0*pff_sum

    f01=(-puf_sum+sqrt(puf_sum**2+4*pff_sum*delta_t*p_aim)) &
            /(2*delta_t*pff_sum)
    f02=(-puf_sum-sqrt(puf_sum**2+4*pff_sum*delta_t*p_aim)) &
            /(2*delta_t*pff_sum)
    if (abs(f01)<abs(f02)) then
        do j=0,tnky
            do k=0,tnkz_s
                do i=0,nkx_s
                    cf1(i,k,j)=f01*cf1(i,k,j)
                    cf3(i,k,j)=f01*cf3(i,k,j)
                end do
            end do
        end do
    else
        do j=0,tnky
            do k=0,tnkz_s
                do i=0,nkx_s
                    cf1(i,k,j)=f02*cf1(i,k,j)
                    cf3(i,k,j)=f02*cf3(i,k,j)
                end do
            end do
        end do
    end if

else if (f_type.eq.3) then
! **** PROPAGATING WAVE FORCING ****
    do j=0,tnky
        do k=0,tnkz_s
            do i=0,nkx_s
                if ((kx_s(i).eq.0) .and. ((kz_s(k).lt.0) .or. &
                        ((kz_s(k).eq.0) .and. (ky(j).lt.0)))) then
                else
                    if ((i+rankz*(nkx_s+1).le.nkx) .and. &
                            (k+ranky*(tnkz_s+1).le.tnkz)) then
                        kappa2=kx2_s(i)+kz2_s(k)
                        k2=kappa2+ky2(j)
    !                if ((k2.le.100.) .and. (ky(j).ne.0)
                        if ((k2.gt.6.25) .and. (k2.lt. 12.25) &
                                .and. (kappa2.ne.0)) then
                            cs1(i,k,j)=cexp(cmplx(0,f_phase(i,k,j) - &
                                        sqrt(ri_tau(1)*kappa2/k2)*time))
    !     &                    *(kappa2**3*k2)**(-0.5)
                            cf1(i,k,j)=cs1(i,k,j)*ky(j)*kx_s(i)/sqrt(k2*kappa2)
                            cf2(i,k,j)=-cs1(i,k,j)*sqrt(kappa2)/sqrt(k2)
                            cf3(i,k,j)=cs1(i,k,j)*ky(j)*kz_s(k)/sqrt(k2*kappa2)
                            csth1(i,k,j)=cs1(i,k,j)*ci/sqrt(ri_tau(1))
                            puf=puf+conjg(cu1(i,k,j))*cf1(i,k,j) &
                                    +conjg(cu2(i,k,j))*cf2(i,k,j) &
                                    +conjg(cu3(i,k,j))*cf3(i,k,j) &
                                    +ri_tau(1)*conjg(cth(i,k,j,1))*csth1(i,k,j)
                            if (first_time) then
                                pff=pff+1.d0 !abs(cs1(i,k,j))**2
                            end if
                        end if
                    end if
                end if
            end do
        end do
    end do

    call mpi_allreduce(puf,puf_sum,1,mpi_double_precision,mpi_sum, &
                        mpi_comm_world,ierror)
    if (first_time) then
        call mpi_allreduce(pff,pff_sum,1,mpi_double_precision, &
                            mpi_sum,mpi_comm_world,ierror)
    end if
    puf_sum=2.d0*puf_sum
    if (first_time) pff_sum=2.d0*pff_sum

    f01=(-puf_sum+sqrt(puf_sum**2+4*pff_sum*delta_t*p_aim)) &
            /(2*delta_t*pff_sum)
    f02=(-puf_sum-sqrt(puf_sum**2+4*pff_sum*delta_t*p_aim)) &
            /(2*delta_t*pff_sum)

    if (abs(f01)<abs(f02)) then
        do j=0,tnky
            do k=0,tnkz_s
                do i=0,nkx_s
                    cf1(i,k,j)=f01*cf1(i,k,j)
                    cf2(i,k,j)=f01*cf2(i,k,j)
                    cf3(i,k,j)=f01*cf3(i,k,j)
                    cfth(i,k,j,1)=cfth(i,k,j,1)+f01*csth1(i,k,j)
                end do
            end do
        end do
    else
        do j=0,tnky
            do k=0,tnkz_s
                do i=0,nkx_s
                    cf1(i,k,j)=f02*cf1(i,k,j)
                    cf2(i,k,j)=f02*cf2(i,k,j)
                    cf3(i,k,j)=f02*cf3(i,k,j)
                    cfth(i,k,j,1)=cfth(i,k,j,1)+f02*csth1(i,k,j)
                end do
            end do
        end do
    end if

end if

return
end
