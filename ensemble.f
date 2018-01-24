C******************************************************************************|
C ensemble.f -> main program for Diablo-EnKF and all required subroutines      
C
C******************************************************************************|

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      program ensemble
	include 'header_ensem'
	integer i,j,k,status(mpi_status_size)
	real*8 simtime
	logical n
	
	call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,rank,ierr)
      call mpi_comm_size(mpi_comm_world,size,ierr)
	if (rank.eq.root .and. NE.ne.size-1) then
	write(6,*)'--------------'
	write(6,*)' Fatal Error:'
	write(6,*)'   Process number inconsistnent with number of'
	write(6,*)'   ensemble members, terminating program.'
	write(6,*)'--------------'
      end if
	if (NE.ne.size-1) STOP 
	
	call initialize

	if (rank.eq.root) then
	  call system('clear;clear',i)
        write(6,*) 
	  write(6,*) '             *******************************'
        write(6,*) '             ****** Welcome to Diablo ******'
	  write(6,*) '             **** Ensemble Kalman Filter ***'
	  write(6,*) '             *******************************'
        write(6,*)
	end if
	
      call mpi_barrier(mpi_comm_world,ierr)
	
      simtime = 0.d0
	time_step = 0
	n = .true.
	do while (simtime.lt.20)
	  time_step = time_step+1
	  simtime=simtime+delta_t
	  if (rank.eq.root) then
	    write(6,'(A,F8.3,I5)') 'Simulation Time: ', simtime,time_step
	  end if
	  do rk_step=1,3
	    call rk_per_1
	  end do
	  
	  if (mod(time_step,meas_freq).eq.0) then
	    if (time_step.gt.10) then
	      if (rank.eq.root) write(6,*) '.....updating ensembles.....'
	      call update
		if (n) then
	        call save_stats(.false.)
		  n = .false.
		else
		  n = .true.
		end if
	    end if
	  end if 
 
	  if (mod(time_step,save_stats_int).eq.0) then
C	    call save_stats(.false.)
	  end if
      end do
	call save_stats(.true.)
	call mpi_finalize(ierr)
	
	if (rank.eq.root) then
        write(6,*)
        write(6,*) '        ****** I''m done ******'
        write(6,*)
	end if
      end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine update
	include 'header_ensem'
	integer i,j,k,v,m,status(mpi_status_size)
	real*8 mu
	
	call average
	sze = NV*NM
	if (rank.ne.root) then
	  do n=1,n_th
	    do j=0,tnky
	      do k=0,tnkz
		  do i=0,nkx
		    CRTH(i,k,j,n) = CTH(i,k,j,n)-CRTH(i,k,j,n)
		    CR1(i,k,j) = CU1(i,k,j)-CR1(i,k,j)
		    CR3(i,k,j) = CU3(i,k,j)-CR3(i,k,j)
	        end do
	      end do
	    end do
	  end do
	  do n=1,3
	    rpos(n) = pos(n)-rpos(n)
	  end do
	  rtheta = theta-rtheta
	  call output(z,veh_pos(:,:,0),.FALSE.,CRTH,CR1,CR3)
	  call distanceR
	  do m=1,sze
	    do v = m,sze
	      Rbuff(v,m) = R(v,m)*z(m,1)*z(v,1) / (NE-1)
		if (v.ne.m) Rbuff(m,v) = Rbuff(v,m)
	    end do
	  end do
	  do m=1,NM
	    do v=1,NV
	      i = NV*(m-1)+v
		Rbuff(i,i) = Rbuff(i,i) + meas_var(m)
	    end do
	  end do
	else
	  Rbuff = 0.
	  Q = 0.
	  gamma = 0.
	  call output(z,veh_pos(:,:,0),.TRUE.,CTH,CU1,CU3)
	end if
	call mpi_allreduce(Rbuff,R,sze*sze,mpi_real8,mpi_sum,
     &                   mpi_comm_world,ierr)

	sze = (3*nx*nz+3)*nv*nm
	if (rank.ne.root) then
	  do n=1,n_th
	    call fft_xzy_to_physical(CRTH(0,0,0,n),RTH(0,0,0,n))
	  end do
	  call fft_xzy_to_physical(CR1,R1)
	  call fft_xzy_to_physical(CR3,R3)
	  call distanceQ
	  do j=1,nv*nm
	    do n=1,n_th
	      do k=0,nzm
	        do i=0,nxm
		    ind = nx*k + i + 1
		    Qbuff(ind,j) = Q(ind,j)*RTH(i,k,0,n)*z(j,1)/(NE-1.)
		    Qbuff(nx*nz+ind,j) = Q(nx*nz+ind,j)*
     &					 R1(i,k,0)*z(j,1)/(NE-1.)
		    Qbuff(2*nx*nz+ind,j) = Q(2*nx*nz+ind,j)*
     &					   R3(i,k,0)*z(j,1)/(NE-1.)
		  end do
	      end do
	    end do
	    if (j.le.nv) then
	      Qbuff(3*nx*nz+1,j) = Q(3*nx*nz+1,j)*rpos(1)*z(j,1) /(NE-1.)
	      Qbuff(3*nx*nz+2,j) = Q(3*nx*nz+2,j)*rpos(3)*z(j,1) /(NE-1.)
	      Qbuff(3*nx*nz+3,j) = Q(3*nx*nz+3,j)*rtheta*z(j,1) /(NE-1.)
	    else
	      Qbuff(3*nx*nz+1,j) = 0.
	      Qbuff(3*nx*nz+2,j) = 0.
	      Qbuff(3*nx*nz+3,j) = 0.
	    end if
	  end do
	end if
	call mpi_allreduce(Qbuff,Q,sze,mpi_real8,mpi_sum,
     &                   mpi_comm_world,ierr)
       
	
      if (rank.ne.root) then
	  mu =0.
	  do m=1,NM
	    do v=1,NV
	      gamma(v,m) = zbqlnor(mu,sqrt(meas_var(m)))
	    end do
	  end do
	end if
	sze = nv*nm
	call mpi_bcast(z,sze,mpi_real8,root,mpi_comm_world,ierr)
	call mpi_allreduce(gamma,Rbuff(:,1),sze,mpi_real8,mpi_sum,
     &			 mpi_comm_world,ierr)
     	call mpi_barrier(mpi_comm_world,ierr)
	
	if (rank.ne.root) then
	  do m=1,NM
	    do v=1,NV
	      gamma(v,m) = gamma(v,m) - Rbuff(NV*(m-1)+v,1)/NE
		z(v,m) = z(v,m) + gamma(v,m)
	    end do
	  end do
	  call output(gamma,veh_pos(:,:,0),.FALSE.,CTH,CU1,CU3)
	  do m=1,NM
	    do v=1,NV
		gamma(v,m) = z(v,m) - gamma(v,m)
	    end do
	  end do
	  call gauss
	end if

	
	return
	end
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine average
	include 'header_ensem'
	integer status(mpi_status_size),i,j,k,n
	   
	sze = (NX+2)*(NZ+2)*(NY+2)*N_TH
	call mpi_reduce(pos,rpos,3,mpi_real8,mpi_sum,
     &                root,mpi_comm_world,ierr)
      call mpi_reduce(theta,rtheta,1,mpi_real8,mpi_sum,
     &                root,mpi_comm_world,ierr)
	call mpi_reduce(TH,RTH,sze,mpi_real8,mpi_sum,
     &                root,mpi_comm_world,ierr)
      call mpi_reduce(U1,R1,sze,mpi_real8,mpi_sum,
     &                root,mpi_comm_world,ierr)
      call mpi_reduce(U3,R3,sze,mpi_real8,mpi_sum,
     &                root,mpi_comm_world,ierr)
	if (rank.eq.root) then
	  do n=1,n_th
	    do j=0,tnky
	      do k=0,tnkz
		  do i=0,nkx
		    CRTH(i,k,j,n) = (CRTH(i,k,j,n)-CTH(i,k,j,n))/NE
		    CR1(i,k,j) = (CR1(i,k,j)-CU1(i,k,j))/NE
		    CR3(i,k,j) = (CR3(i,k,j)-CU3(i,k,j))/NE
	        end do
	      end do
	    end do
	  end do
	  do n=1,3
	    rpos(n) = (rpos(n)-pos(n))/NE
	  end do
	  rtheta = (rtheta-theta)/NE
	end if
	
	call mpi_bcast(rpos,3,mpi_real8,root,mpi_comm_world,ierr)
	call mpi_bcast(rtheta,1,mpi_real8,root,mpi_comm_world,ierr)
	call mpi_bcast(RTH,sze,mpi_real8,
     &		        root,mpi_comm_world,ierr)
      call mpi_bcast(R1,sze,mpi_real8,
     &		        root,mpi_comm_world,ierr)
      call mpi_bcast(R3,sze,mpi_real8,
     &		        root,mpi_comm_world,ierr)
     
	return
	end
	

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine output(YZ,veh,noise,TCTH,TCU1,TCU3)
	include 'header_ensem'
	real*8 veh(1:3,1:NV), zbqlnor, mu, YZ(1:NV,1:NM)
	complex*16 TCTH(0:NX/2,0:NZ+1,0:NY+1,1:N_TH),
     &	     TCU1(0:NX/2,0:NZ+1,0:NY+1),
     &	     TCU3(0:NX/2,0:NZ+1,0:NY+1)
	logical noise
	integer n,m,v
		
	if (num_per_dir.eq.3) then
	
	  YZ = 0.
	  do m=1,NM
	    select case(m)
	      case(1)
		  do v=1,NV
		    do n=1,N_TH
	          call meas_per(TCTH(:,:,:,n),veh(:,v),YZ(v,m))
		    end do
	        end do
		case(2)
		  do v=1,NV
	          call meas_per(TCU1,veh(:,v),YZ(v,m))
	        end do
	      case(3)
		  do v=1,NV
	          call meas_per(TCU3,veh(:,v),YZ(v,m))
	        end do
		case default
		  write(6,*) 'Measurement type not supported.'
	    end select
	  end do

	else
	  write(6,*) 'Non-periodic cases not supported'
	end if
	
	if (noise) then
	  mu =0.
	  do m=1,NM
	    do v=1,NV
	      YZ(v,m) = YZ(v,m) + zbqlnor(mu,sqrt(meas_var(m)))
	    end do
	  end do
	end if
	
	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine meas_per(CU,veh,zi)
	include 'header_ensem'
	integer i,k,j
      complex*16 CU(0:NX/2,0:NZ+1,0:NY+1)
	real*8 veh(3),zi
	
	zi = 0.
	do j=0,tnky
	  do k=0,tnkz
	    do i=0,nkx
 		zi = zi + CU(i,k,j) *
     &	exp(-meas_ave**2*(KX2(i)+KY2(j)+KZ2(k))/2.d0) *
     &	exp(-ci*(KX(i)*(LX-veh(1))
     &             + KY(j)*(LY-veh(2))
     &             + KZ(k)*(LZ-veh(3))))
		if (i.gt.0) then
		zi = zi + conjg(CU(i,k,j)) *
     &	exp(-meas_ave**2*(KX2(i)+KY2(j)+KZ2(k))/2.d0) *
     &	exp(ci*(KX(i)*(LX-veh(1))
     &            + KY(j)*(LY-veh(2))
     &            + KZ(k)*(LZ-veh(3))))
		end if
	    end do
	  end do
	end do

	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine gauss
	include 'header_ensem'
	integer n,i,j,k,ind
	
	sze = NV*NM
	do j=1,sze-1
	  do i=j+1,sze
	    R(i,j) = -R(i,j) / R(j,j)
	  end do
	  do i=j+1,sze
	    do k=j+1,sze
	      R(k,i) = R(k,i) + R(k,j)*R(j,i)
	    end do
	  end do
	  do i=j+1,sze
	    gamma(i,1) = gamma(i,1) + R(i,j) * gamma(j,1)
	  end do
	end do
	gamma(sze,1) = gamma(sze,1) / R(sze,sze)
	do i=sze-1,1,-1
	  do j=i+1,sze
	    gamma(i,1) = gamma(i,1) - R(i,j)*gamma(j,1)
	  end do
	  gamma(i,1) = gamma(i,1) / R(i,i)
	end do
	
	do n=1,n_th
	  call fft_xzy_to_physical(CTH(0,0,0,n),TH(0,0,0,n))
	end do
	call fft_xzy_to_physical(CU1,U1)
	call fft_xzy_to_physical(CU3,U3)
	do j=1,nv*nm
	  do n=1,n_th
	    do k=0,nzm
	      do i=0,nxm
	        ind = nx*k + i + 1
		  TH(i,k,0,n) = TH(i,k,0,n) + Q(ind,j)*gamma(j,1)
		  U1(i,k,0) = U1(i,k,0) + Q(nx*nz+ind,j)*gamma(j,1)
		  U3(i,k,0) = U3(i,k,0) + Q(2*nx*nz+ind,j)*gamma(j,1)
		end do
	    end do
	  end do
	  pos(1) = pos(1) + Q(3*nx*nz+1,j)*gamma(j,1)
	  pos(3) = pos(3) + Q(3*nx*nz+2,j)*gamma(j,1)
	  theta  = theta  + Q(3*nx*nz+3,j)*gamma(j,1)
	end do
	
	do n=1,n_th
	  call fft_xzy_to_fourier(TH(0,0,0,n),CTH(0,0,0,n))
	end do
	call fft_xzy_to_fourier(U1,CU1)
	call fft_xzy_to_fourier(U3,CU3)

	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine distanceR
	include 'header_ensem'
	integer m,v,i,j,ind1,ind2
	real*8 xd,zd,rho
	
	do m=1,NM
	  do v=1,NV
	    do i=1,NM
	      do j=1,NV
		  ind1 = (i-1)*NV+j
		  ind2 = (m-1)*NV+v
	        xd = veh_pos(1,v,0)-veh_pos(1,j,0)
		  zd = veh_pos(3,v,0)-veh_pos(3,j,0)
		  xd = sqrt(xd**2+zd**2)
		  R(ind1,ind2) = rho(xd)
	      end do
	    end do
	  end do
	end do
	
	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine distanceQ
	include 'header_ensem'
	integer m,v,i,k,ind1,ind2
	real*8 xd,zd,rho
	
	do m=1,NM
	  do v=1,NV
	    do i=0,nxm
	      do k=0,nzm
		  ind1 = nx*k + i + 1
		  ind2 = (m-1)*NV + v
		  xd = veh_pos(1,v,0) - gx(i)
		  zd = veh_pos(3,v,0) - gz(k)
		  xd = sqrt(xd**2+zd**2)
		  Q(ind1,ind2) = rho(xd)
		  Q(nx*nz+ind1,ind2) = rho(xd)
		  Q(2*nx*nz+ind1,ind2) = rho(xd)
		end do
	    end do
	    xd = veh_pos(1,v,0)-pos(1)
	    zd = veh_pos(3,v,0)-pos(3)
	    xd = sqrt(xd**2+zd**2)
	    Q(3*nx*nz+1,ind2) = rho(xd*2.)
	    Q(3*nx*nz+2,ind2) = rho(xd*2.)
	    Q(3*nx*nz+3,ind2) = rho(xd*2.)
	  end do
	end do
	
	
	
	return
	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	function rho(dd)
	real*8 rho,dd,rad
	
	rad = 0.5
	rho = exp(-((2.*dd/rad)**2))

	end
	
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
	subroutine variance
	include 'header_ensem'
	integer n,i,k
	
	sze = (NX+2)*(NZ+2)*(NY+2)*N_TH
	if (rank.ne.root) then
	  do n=1,n_th
	    call fft_xzy_to_physical(CTH(0,0,0,n),TH(0,0,0,n))
	    call fft_xzy_to_physical(CRTH(0,0,0,n),RTH(0,0,0,n))
	  end do
	  call fft_xzy_to_physical(CU1,U1)
	  call fft_xzy_to_physical(CU3,U3)
	  call fft_xzy_to_physical(CR1,R1)
	  call fft_xzy_to_physical(CR3,R3)
	  do n=1,n_th
	    do k=0,nzm
	      do i=0,nxm
		  FTH(i,k,0,n) = ((TH(i,k,0,n)-RTH(i,k,0,n))**2)/(NE-1)
		  F1(i,k,0) = ((U1(i,k,0)-R1(i,k,0))**2)/(NE-1)
		  F3(i,k,0) = ((U3(i,k,0)-R3(i,k,0))**2)/(NE-1)
	      end do
	    end do
	  end do
	  do n=1,n_th
	    call fft_xzy_to_fourier(TH(0,0,0,n),CTH(0,0,0,n))
	  end do
	  call fft_xzy_to_fourier(U1,CU1)
	  call fft_xzy_to_fourier(U3,CU3)
	else
	  FTH = 0.
	  F1 = 0.
	  F3 = 0.
	end if
	call mpi_reduce(FTH(:,:,:,1),F2,sze,mpi_real8,mpi_sum,
     &                root,mpi_comm_world,ierr)
      if (rank.eq.root) then
	  do n=1,n_th
	    do k=0,nzm
	      do i=0,nxm
		  FTH(i,k,0,n) = F2(i,k,0)
	      end do
	    end do
	  end do
	end if
	call mpi_reduce(F1,F2,sze,mpi_real8,mpi_sum,
     &                root,mpi_comm_world,ierr)
      if (rank.eq.root) then
	  do k=0,nzm
	    do i=0,nxm
	      F1(i,k,0) = F2(i,k,0)
	    end do
	  end do
	end if
	call mpi_reduce(F3,F2,sze,mpi_real8,mpi_sum,
     &                root,mpi_comm_world,ierr)
      if (rank.eq.root) then
	  do k=0,nzm
	    do i=0,nxm
	      F3(i,k,0) = F2(i,k,0)
	    end do
	  end do
	end if
	
	
	
	return
	end

