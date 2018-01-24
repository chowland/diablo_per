C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE NETCDF_OPEN_VIS
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE 'netcdf.inc'

      integer dim_3d_1(3),dim_3d_2(3),dim_3d_3(3),dim_3d_4(3)
      integer x_dimid, y_dimid, z_dimid, x_varid, y_varid, z_varid
      integer xf_dimid,yf_dimid,zf_dimid,xfvarid,yf_varid,zf_varid
      integer n_th_dimid,n_th_varid
      integer dim_th(4)
      integer retval
      integer N_TH_ARRAY(1:N_TH)

      integer i,j,n

! Create a new netCDF file, overwrite if it already exists     
      retval = nf_create('vis.nc', NF_CLOBBER, ncid_vis)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     Define the dimensions. NetCDF will hand back an ID for each.
      retval = nf_def_dim(ncid_vis, "x", NX, x_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid_vis, "y", NY, y_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid_vis, "z", NZ, z_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid_vis, "xf", NX, xf_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid_vis, "yf", NY, yf_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid_vis, "zf", NZ, zf_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid_vis, "n_th", N_TH, n_th_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)

! Define the coordinate variables in the netCDF file
      retval=nf_def_var(ncid_vis,"x",NF_DOUBLE,1,x_dimid,x_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_vis,"y",NF_DOUBLE,1,y_dimid,y_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_vis,"z",NF_DOUBLE,1,z_dimid,z_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_vis,"xf",NF_DOUBLE,1,xf_dimid,xf_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_vis,"yf",NF_DOUBLE,1,yf_dimid,yf_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_vis,"zf",NF_DOUBLE,1,zf_dimid,zf_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_vis,"n_th",NF_INT,1,n_th_dimid,n_th_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     The dimids array is used to pass the IDs of the dimensions of
C     the variables. Note that in fortran arrays are stored in
C     column-major format.

      IF (NUM_PER_DIR.eq.3) THEN
        dim_3d_1(1) = x_dimid
        dim_3d_1(2) = y_dimid
        dim_3d_1(3) = z_dimid

        dim_3d_2(1) = x_dimid
        dim_3d_2(2) = y_dimid
        dim_3d_2(3) = z_dimid

        dim_3d_3(1) = x_dimid
        dim_3d_3(2) = y_dimid
        dim_3d_3(3) = z_dimid

        dim_3d_4(1) = x_dimid
        dim_3d_4(2) = y_dimid
        dim_3d_4(3) = z_dimid

        dim_th(1) = x_dimid
        dim_th(2) = y_dimid
        dim_th(3) = z_dimid
        dim_th(4) = n_th_dimid
      ELSE IF (NUM_PER_DIR.eq.2) THEN
        dim_3d_1(1) = x_dimid
        dim_3d_1(2) = yf_dimid
        dim_3d_1(3) = z_dimid
      
        dim_3d_2(1) = x_dimid
        dim_3d_2(2) = y_dimid
        dim_3d_2(3) = z_dimid

        dim_3d_3(1) = x_dimid
        dim_3d_3(2) = yf_dimid
        dim_3d_3(3) = z_dimid

        dim_3d_4(1) = x_dimid
        dim_3d_4(2) = yf_dimid
        dim_3d_4(3) = z_dimid

        dim_th(1) = x_dimid
        dim_th(2) = y_dimid
        dim_th(3) = z_dimid
        dim_th(4) = n_th_dimid
      ELSE IF (NUM_PER_DIR.eq.1) THEN
        dim_3d_1(1) = x_dimid
        dim_3d_1(2) = yf_dimid
        dim_3d_1(3) = zf_dimid
        
        dim_3d_2(1) = x_dimid
        dim_3d_2(2) = y_dimid
        dim_3d_2(3) = zf_dimid

        dim_3d_3(1) = x_dimid
        dim_3d_3(2) = yf_dimid
        dim_3d_3(3) = z_dimid

        dim_3d_4(1) = x_dimid
        dim_3d_4(2) = yf_dimid
        dim_3d_4(3) = zf_dimid

        dim_th(1) = x_dimid
        dim_th(2) = yf_dimid
        dim_th(3) = zf_dimid
        dim_th(4) = n_th_dimid
      ELSE 
! Else, we must be in the Cavity flow case
        dim_3d_1(1) = x_dimid
        dim_3d_1(2) = yf_dimid
        dim_3d_1(3) = zf_dimid

        dim_3d_2(1) = xf_dimid
        dim_3d_2(2) = y_dimid
        dim_3d_2(3) = zf_dimid

        dim_3d_3(1) = xf_dimid
        dim_3d_3(2) = yf_dimid
        dim_3d_3(3) = z_dimid

        dim_3d_4(1) = xf_dimid
        dim_3d_4(2) = yf_dimid
        dim_3d_4(3) = zf_dimid

        dim_th(1) = xf_dimid
        dim_th(2) = yf_dimid
        dim_th(3) = zf_dimid
        dim_th(4) = n_th_dimid
      END IF

! Define the variable data in our netCDF file
      retval=nf_def_var(ncid_vis,"U",NF_DOUBLE,3,dim_3d_1
     &        ,u_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_vis,"V",NF_DOUBLE,3,dim_3d_2
     &        ,v_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_vis,"W",NF_DOUBLE,3,dim_3d_3
     &        ,w_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_vis,"P",NF_DOUBLE,3,dim_3d_4
     &        ,p_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_vis,"TH",NF_DOUBLE,4,dim_th
     &        ,th_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     End define mode. This tells netCDF we are done defining metadata.
      retval = nf_enddef(ncid_vis)
      if (retval .ne. nf_noerr) call handle_err(retval)
 
! The base grid points will be indexed from 0:N-1 if the direction is
! periodic, and from 1:N if the direction is wall-bounded:
      IF (NUM_PER_DIR.eq.3) THEN
        retval = nf_put_var_double(ncid_vis,x_varid,GX(0:NX-1))
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid_vis,y_varid,GY(0:NY-1))
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid_vis,z_varid,GZ(0:NZ-1))
        if (retval .ne. nf_noerr) call handle_err(retval)
      ELSE IF (NUM_PER_DIR.eq.2) THEN
        retval = nf_put_var_double(ncid_vis,x_varid,GX(0:NX-1))
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid_vis,y_varid,GY(1:NY))
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid_vis,z_varid,GZ(0:NZ-1))
        if (retval .ne. nf_noerr) call handle_err(retval)
      ELSE IF (NUM_PER_DIR.eq.1) THEN
        retval = nf_put_var_double(ncid_vis,x_varid,GX(0:NX-1))
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid_vis,y_varid,GY(1:NY))
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid_vis,z_varid,GZ(1:NZ))
        if (retval .ne. nf_noerr) call handle_err(retval)
      ELSE
        retval = nf_put_var_double(ncid_vis,x_varid,GX(1:NX))
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid_vis,y_varid,GY(1:NY))
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid_vis,z_varid,GZ(1:NZ))
        if (retval .ne. nf_noerr) call handle_err(retval)
      END IF

! The fractional grid is always indexed from 1:N
        retval = nf_put_var_double(ncid_vis,xf_varid,GXF(1:NX))
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid_vis,yf_varid,GYF(1:NY))
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_double(ncid_vis,zf_varid,GZF(1:NZ))
        if (retval .ne. nf_noerr) call handle_err(retval)

! Create an array of N_TH indexes
      do n=1,N_TH
        N_TH_ARRAY(n)=n
      end do
      retval = nf_put_var_int(ncid_vis,n_th_varid,N_TH_ARRAY(1:N_TH))
      if (retval.ne.nf_noerr) call handle_err(retval)

! Define the start and count arrays
      nc_count(1)=NY
      nc_count(2)=1
      nc_start(1)=1 
      nc_index=0
 
      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE NETCDF_OPEN_STATS_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE 'netcdf.inc'

      integer dim_2d(2)
      integer dim_2d_t(3), dim_1d_t(2)
      integer x_dimid, y_dimid, z_dimid, x_varid, y_varid, z_varid
      integer t_dimid
      integer retval

      integer i,j


! Create a new netCDF file, overwrite if it already exists     
      retval = nf_create('mean.nc', NF_CLOBBER, ncid_stats)
      if (retval .ne. nf_noerr) call handle_err(retval)


!     Define the dimensions. NetCDF will hand back an ID for each.
      retval = nf_def_dim(ncid_stats, "x", NX, x_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid_stats, "y", NY, y_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid_stats, "z", NZ, z_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid_stats, "t", NF_UNLIMITED, t_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)

! Define the coordinate variables in the netCDF file
      retval=nf_def_var(ncid_stats,"x",NF_DOUBLE,1,x_dimid,x_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"y",NF_DOUBLE,1,y_dimid,y_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"z",NF_DOUBLE,1,z_dimid,z_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"t",NF_DOUBLE,1,t_dimid,t_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     The dimids array is used to pass the IDs of the dimensions of
C     the variables. Note that in fortran arrays are stored in
C     column-major format.

      dim_2d(1) = x_dimid
      dim_2d(2) = y_dimid

      dim_2d_t(1) = x_dimid
      dim_2d_t(2) = y_dimid
      dim_2d_t(3) = t_dimid

      dim_1d_t(1) = y_dimid
      dim_1d_t(2) = t_dimid

! Define the variable data in our netCDF file
      retval=nf_def_var(ncid_stats,"U_BAR",NF_DOUBLE,2,dim_1d_t
     &        ,ubar_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"V_BAR",NF_DOUBLE,2,dim_1d_t
     &        ,vbar_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"W_BAR",NF_DOUBLE,2,dim_1d_t
     &        ,wbar_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"U_RMS",NF_DOUBLE,2,dim_1d_t
     &        ,urms_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"V_RMS",NF_DOUBLE,2,dim_1d_t
     &        ,vrms_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"W_RMS",NF_DOUBLE,2,dim_1d_t
     &        ,wrms_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"uv",NF_DOUBLE,2,dim_1d_t
     &        ,uv_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"wv",NF_DOUBLE,2,dim_1d_t
     &        ,wv_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval=nf_def_var(ncid_stats,"uw",NF_DOUBLE,2,dim_1d_t
     &        ,uw_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     End define mode. This tells netCDF we are done defining metadata.
      retval = nf_enddef(ncid_stats)
      if (retval .ne. nf_noerr) call handle_err(retval)
    
      retval = nf_put_var_double(ncid_stats,x_varid,GX(0:NX-1))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_double(ncid_stats,z_varid,GZ(0:NZ-1))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_double(ncid_stats,y_varid,GYF(1:NY))
      if (retval .ne. nf_noerr) call handle_err(retval)

! Define the start and count arrays
      nc_count(1)=NY
      nc_count(2)=1
      nc_start(1)=1 
      nc_index=0
 
      RETURN
      END

      subroutine NETCDF_WRITE_STATS_CHAN
      include 'header'
      include 'netcdf.inc'
      integer retval

      nc_index=nc_index+1
      nc_start(2)=nc_index
C Write the data to the netCDF file
      retval = nf_put_vara_double(ncid_stats,urms_varid
     &             ,nc_start,nc_count,urms(1:NY))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,vrms_varid
     &             ,nc_start,nc_count,vrms(1:NY))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,wrms_varid
     &             ,nc_start,nc_count,wrms(1:NY))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,uv_varid
     &             ,nc_start,nc_count,uv(1:NY))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,uw_varid
     &             ,nc_start,nc_count,uw(1:NY))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,wv_varid
     &             ,nc_start,nc_count,wv(1:NY))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,ubar_varid
     &             ,nc_start,nc_count,dble(CR1(0,0,1:NY)))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,vbar_varid
     &             ,nc_start,nc_count,dble(CR2(0,0,1:NY)))
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_vara_double(ncid_stats,wbar_varid
     &             ,nc_start,nc_count,dble(CR3(0,0,1:NY)))
      if (retval .ne. nf_noerr) call handle_err(retval)

C Write the current time to the netCDF file
      retval = nf_put_var_double(ncid_stats,t_varid,TIME)
      if (retval .ne. nf_noerr) call handle_err(retval)

      return
      end

      subroutine NETCDF_WRITE_VIS
      include 'header'
      include 'netcdf.inc'
      integer retval
      integer i,j,k
      real*8 tempu(NX,NY,NZ)
      real*8 temp_th(NX,NY,NZ,1:N_TH)

      integer ii1,ii2,jj1,jj2,kk1,kk2
  
! Set the bounds for loops 
      IF (NUM_PER_DIR.eq.0) THEN
        ii1=1
        ii2=NX
        jj1=1
        jj2=NY
        kk1=1
        kk2=NZ
      ELSE IF (NUM_PER_DIR.eq.1) THEN
        ii1=0
        ii2=NX-1
        jj1=1
        jj2=NY
        kk1=1
        kk2=NZ
      ELSE IF (NUM_PER_DIR.eq.2) THEN
        ii1=0
        ii2=NX-1
        jj1=1
        jj2=NY
        kk1=0
        kk2=NZ-1
      ELSE IF (NUM_PER_DIR.eq.3) THEN
        ii1=0
        ii2=NX-1
        jj1=0
        jj2=NY-1
        kk1=0
        kk2=NZ-1
      END IF 

      IF (NUM_PER_DIR.eq.2) THEN
      call FFT_XZ_TO_PHYSICAL(CU1,U1,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU2,U2,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CU3,U3,0,NY+1)
      call FFT_XZ_TO_PHYSICAL(CP,P,0,NY+1)
      DO N=1,N_TH
        CALL FFT_XZ_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N),0,NY+1)
      END DO
      ELSE IF (NUM_PER_DIR.eq.3) THEN
      call FFT_XZY_TO_PHYSICAL(CU1,U1)
      call FFT_XZY_TO_PHYSICAL(CU2,U2)
      call FFT_XZY_TO_PHYSICAL(CU3,U3)
      call FFT_XZY_TO_PHYSICAL(CP,P)
      DO N=1,N_TH
        CALL FFT_XZY_TO_PHYSICAL(CTH(0,0,0,N),TH(0,0,0,N))
      END DO
      END IF

      do i=ii1,ii2
      do k=kk1,kk2
      do j=jj1,jj2
        tempu(i+1-ii1,j+1-jj1,k+1-kk1)=U1(i,k,j)
      end do
      end do
      end do
C Write the 3d velocity field:
      retval=nf_put_var_double(ncid_vis,u_varid
     &         ,tempu)
      if (retval .ne. nf_noerr) call handle_err(retval)

      do i=ii1,ii2
      do k=kk1,kk2
      do j=jj1,jj2
        tempu(i+1-ii1,j+1-jj1,k+1-kk1)=U2(i,k,j)
      end do
      end do
      end do
C Write the 3d velocity field:
      retval=nf_put_var_double(ncid_vis,v_varid
     &         ,tempu)
      if (retval .ne. nf_noerr) call handle_err(retval)

      do i=ii1,ii2
      do k=kk1,kk2
      do j=jj1,jj2
        tempu(i+1-ii1,j+1-jj1,k+1-kk1)=U3(i,k,j)
      end do
      end do
      end do
C Write the 3d velocity field:
      retval=nf_put_var_double(ncid_vis,w_varid
     &         ,tempu)
      if (retval .ne. nf_noerr) call handle_err(retval)

      do i=ii1,ii2
      do k=kk1,kk2
      do j=jj1,jj2
        tempu(i+1-ii1,j+1-jj1,k+1-kk1)=P(i,k,j)
      end do
      end do
      end do
C Write the 3d pressure field:
      retval=nf_put_var_double(ncid_vis,p_varid
     &         ,tempu)
      if (retval .ne. nf_noerr) call handle_err(retval)

      do i=ii1,ii2
      do k=kk1,kk2
      do j=jj1,jj2
      do n=1,N_TH
        temp_th(i+1-ii1,j+1-jj1,k+1-kk1,n)=TH(i,k,j,n)
      end do
      end do
      end do
      end do
C Write the 3d scalar field:
      retval=nf_put_var_double(ncid_vis,th_varid
     &         ,temp_th)
      if (retval .ne. nf_noerr) call handle_err(retval)


      return 
      end


      subroutine NETCDF_CLOSE_VIS
      include 'header'
      include 'netcdf.inc'
      integer retval

      IF (NUM_PER_DIR.eq.2) THEN
      call FFT_XZY_TO_FOURIER(U1,CU1)
      call FFT_XZY_TO_FOURIER(U2,CU2)
      call FFT_XZY_TO_FOURIER(U3,CU3)
      call FFT_XZY_TO_FOURIER(P,CP)
      DO N=1,N_TH
        CALL FFT_XZY_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N))
      END DO
      ELSE IF (NUM_PER_DIR.eq.3) THEN
      call FFT_XZY_TO_FOURIER(U1,CU1)
      call FFT_XZY_TO_FOURIER(U2,CU2)
      call FFT_XZY_TO_FOURIER(U3,CU3)
      call FFT_XZY_TO_FOURIER(P,CP)
      DO N=1,N_TH
        CALL FFT_XZY_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N))
      END DO
      END IF

C     Close the file. This frees up any internal netCDF resources
C     associated with the file, and flushes any buffers.
      retval = nf_close(ncid_vis)
      if (retval .ne. nf_noerr) call handle_err(retval)
    
      return
      end

      subroutine NETCDF_CLOSE_STATS_CHAN
      include 'header'
      include 'netcdf.inc'
      integer retval

C     Close the file. This frees up any internal netCDF resources
C     associated with the file, and flushes any buffers.
      retval = nf_close(ncid_stats)
      if (retval .ne. nf_noerr) call handle_err(retval)
    
      return
      end
     

     


      subroutine handle_err(errcode)
      implicit none
      include 'netcdf.inc'
      integer errcode
     
      print *, 'Error: ', nf_strerror(errcode)
      stop 2
      end


