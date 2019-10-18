import numpy as np
import h5py

# Read the dimensionless parameters from the input file
def read_input(rundir):
    fname=rundir+'input.dat'
    f=open(fname,'r')
    for i, line in enumerate(f):
        if i==6:
            NU=float(line.split()[0])
        elif i==19:
            Ri_t=float(line.split()[0])
            Pr=float(line.split()[1])
    f.close()
    Re=1/NU
    return (Re, Ri_t, Pr)

# Read volume averaged quantities from stats.h5
def read_stat(rundir):
    fname=rundir+'stats.h5'
    f=h5py.File(fname,'r')
    nk=f['/U1rms'].attrs.__getitem__('Samples')
    tii=np.zeros(nk)
    U1rms, U2rms, U3rms = np.zeros(nk), np.zeros(nk), np.zeros(nk)
    THrms, THflux = np.zeros(nk), np.zeros(nk)
    epsilon, chi = np.zeros(nk), np.zeros(nk)
    for i in range(nk):
        num=format(i,"04")
        tii[i]=f['/U1rms/'+num].attrs.__getitem__('Time')
        U1rms[i]=f['/U1rms/'+num][()]
        U2rms[i]=f['/U2rms/'+num][()]
        U3rms[i]=f['/U3rms/'+num][()]
        THrms[i]=f['/THrms/'+num][()]
        THflux[i]=f['/THflux/'+num][()]
        epsilon[i]=f['/epsilon/'+num][()]
        chi[i]=f['/chi/'+num][()]
    return (tii, U1rms, U2rms, U3rms, THrms, THflux, epsilon, chi)

# Read horizontally-averaged profiles from mean.h5
def read_mean(rundir):
    fname=rundir+'mean.h5'
    f=h5py.File(fname,'r')
    nk=f['/U1me'].attrs.__getitem__('Samples')
    ny=f['/U1me/0000'].size
    U1me, U3me, THme = np.zeros((nk,ny)), np.zeros((nk,ny)), np.zeros((nk,ny))
    epsilon, chi, THflux = np.zeros((nk,ny)), np.zeros((nk,ny)), np.zeros((nk,ny))
    U1U2, U3U2 = np.zeros((nk,ny)), np.zeros((nk,ny))
    U1U1, U2U2, U3U3, THTH = np.zeros((nk,ny)), np.zeros((nk,ny)), np.zeros((nk,ny)), np.zeros((nk,ny))
    for i in range(nk):
        num=format(i,"04")
        U1me[i]=f['/U1me/'+num][()]
        U3me[i]=f['/U3me/'+num][()]
        THme[i]=f['/THme/'+num][()]
        epsilon[i]=f['/epsilon/'+num][()]
        chi[i]=f['/chi/'+num][()]
        THflux[i]=f['/THflux/'+num][()]
        U1U2[i]=f['/U1U2/'+num][()]
        U3U2[i]=f['/U3U2/'+num][()]
        U1U1[i]=f['/U1rms/'+num][()]
        U2U2[i]=f['/U2rms/'+num][()]
        U3U3[i]=f['/U3rms/'+num][()]
        THTH[i]=f['/THrms/'+num][()]
    return (U1me, U3me, THme, epsilon, chi, THflux, U1U2, U3U2, U1U1, U2U2, U3U3, THTH, nk, ny)

# Read vertical wavenumber spectra from spectra.h5
def read_spectra(rundir):
    fname=rundir+'spectra.h5'
    f=h5py.File(fname,'r')
    nk=f['/U1'].attrs.__getitem__('Samples')
    nky=(f['/U1/0000'].size-1)//2
    U1, U2 = np.zeros((nk,2*nky+1)), np.zeros((nk,2*nky+1))
    U3, TH = np.zeros((nk,2*nky+1)), np.zeros((nk,2*nky+1))
    for i in range(nk):
        num = format(i,"04")
        U1[i]=f['/U1/'+num][()]
        U2[i]=f['/U2/'+num][()]
        U3[i]=f['/U3/'+num][()]
        TH[i]=f['/TH1/'+num][()]
    E1, E2 = np.zeros((nk,nky+1)), np.zeros((nk,nky+1))
    E3, ETH = np.zeros((nk,nky+1)), np.zeros((nk,nky+1))
    E1[:,0], E2[:,0] = U1[:,0], U2[:,0]
    E3[:,0], ETH[:,0] = U3[:,0], TH[:,0]
    for i in range(nky):
        E1[:,i]=U1[:,i]+U1[:,-i]
        E2[:,i]=U2[:,i]+U2[:,-i]
        E3[:,i]=U3[:,i]+U3[:,-i]
        ETH[:,i]=TH[:,i]+TH[:,-i]
    return (E1, E2, E3, ETH, nky)

# Read individual vertical profile from movie_xy.h5 (x-position <-> idx)
def read_profiles(rundir,idx):
    fname=rundir+'movie_xy.h5'
    f=h5py.File(fname,'r')
    nk=f.attrs.__getitem__('Samples')
    ny=f['/0000/U1'].shape[1]
    U1p, U2p, U3p, THp = np.zeros((nk,ny)), np.zeros((nk,ny)), np.zeros((nk,ny)), np.zeros((nk,ny))
    for i in range(nk):
        num = format(i,"04")
        G=f[num+'/U1/'][()]
        U1p[i]=G[:,idx].T
        G=f[num+'/U2/'][()]
        U2p[i]=G[:,idx].T
        G=f[num+'/U3/'][()]
        U3p[i]=G[:,idx].T
        G=f[num+'/TH1/'][()]
        THp[i]=G[:,idx].T
    return (U1p, U2p, U3p, THp)

# import numpy.fft as ft
# import pyfftw



# def out2chi(rundir,outnum):
#     if outnum==0:
#         fname=rundir+'start.h5'
#     elif outnum<10:
#         fname=rundir+'restart_files/out0'+str(outnum)+'.h5'
#     else:
#         fname=rundir+'restart_files/out'+str(outnum)+'.h5'
#     f=h5py.File(fname,'r')
#     S1=f['/TH1/'][()]
#     (NX,NY,NZ)=S1.shape
    
#     Sme=np.mean(S1,axis=(0,2))
#     for i in range(NY):
#         S1[:,i,:]=S1[:,i,:]-Sme[i]
    
#     chi=np.zeros((NX,NY,NZ))
    
#     print('Loaded TH variable')
    
#     A=pyfftw.empty_aligned((NX,NY,NZ),dtype='float64')
#     B1=pyfftw.empty_aligned((NX//2+1,NY,NZ),dtype='complex128')
#     B2=pyfftw.empty_aligned((NX,NY//2+1,NZ),dtype='complex128')
#     B3=pyfftw.empty_aligned((NX,NY,NZ//2+1),dtype='complex128')
    
#     ft1=pyfftw.FFTW(A, B1, axes=(0,), direction='FFTW_FORWARD',threads=16)
#     ift1=pyfftw.FFTW(B1, A, axes=(0,), direction='FFTW_BACKWARD',threads=16)
#     ft2=pyfftw.FFTW(A, B2, axes=(1,), direction='FFTW_FORWARD',threads=16)
#     ift2=pyfftw.FFTW(B2, A, axes=(1,), direction='FFTW_BACKWARD',threads=16)
#     ft3=pyfftw.FFTW(A, B3, axes=(2,), direction='FFTW_FORWARD',threads=16)
#     ift3=pyfftw.FFTW(B3, A, axes=(2,), direction='FFTW_BACKWARD',threads=16)
    
#     CS1=ft1(S1)
#     CIKX=1j*np.arange(NX//2+1).reshape(NX//2+1,1,1)
#     CS1=CIKX*CS1
#     chi=chi+ift1(CS1)
    
#     print('Computed 1st derivative')
    
#     CS1=ft2(S1)
#     CIKY=1j*np.arange(NY//2+1).reshape(1,NY//2+1,1)
#     CS1=CIKY*CS1
#     chi=chi+ift2(CS1)
    
#     print('Computed 2nd derivative')
    
#     CS1=ft3(S1)
#     CIKZ=1j*np.arange(NZ//2+1).reshape(1,1,NZ//2+1)
#     CS1=CIKZ*CS1
#     chi=chi+ift3(CS1)
    
#     print('Computed 3rd derivative')
    
#     (Ri_t, Re, Pr) = (1, 1e4, 1)
#     chi=Ri_t/Re/Pr*chi
    
#     return chi

# # THIS FUNCTION DOESN'T WORK AT THE MOMENT
# # Compute spectra from full output files
# def out2spec(rundir,outnum,Ri_t):
#     if outnum==0:
#         fname=rundir+'start.h5'
#     elif outnum<10:
#         fname=rundir+'restart_files/out0'+str(outnum)+'.h5'
#     else:
#         fname=rundir+'out'+str(outnum)+'.h5'
#     comps=('U1','U2','U3','TH1')
#     f=h5py.File(fname,'r')
#     S1=f['/U2/'][()]
#     Kh2=Ri_t/np.mean(S1**2)
#     (NX,NY,NZ)=S1.shape
#     A=pyfftw.empty_aligned((NX,NY,NZ),dtype='float64')
#     B=pyfftw.empty_aligned((NX,NY,NZ//2+1),dtype='complex128')
#     fft_forw=pyfftw.FFTW(A, B, axes=(0,1,2), direction='FFTW_FORWARD',threads=16)
#     fft_back=pyfftw.FFTW(B, A, axes=(0,1,2), direction='FFTW_BACKWARD',threads=16)
#     spectra={}
#     for var in comps:
#         spectra[var]={}
#         S1=f['/'+var+'/'][()]
#         print('Loaded variable '+var)
        
#         CS1=fft_forw(S1)/(NX*NY*NZ)
#         del S1
#         # UP TO HERE THOUGHT-WISE
#         NKY=int(NY/3)
#         KY=np.arange(NKY+1)
        
#         # Calculate the 2-sided energy spectrum
#         CS1=0.5*CS1*np.conj(CS1)
#         # Dealias the high wavenumbers
#         CS1[NKY+1:-NKY,:,:]=0
#         CS1[:,NKY+1:-NKY,:]=0
#         CS1[:,:,NKY+1:-NKY]=0
        
#         # Calculate the 1-sided energy spectrum
#         CS2=np.zeros((NX,NKY+1,NZ))
#         CS2[:,0,:]=CS1[:,0,:]
#         for j in range(1,NKY+1):
#             CS2[:,j,:]=CS1[:,j,:]+CS1[:,-j,:]
        
#         if var=="TH1":
#             CS2=Ri_t*CS2
        
#         KX, KZ = ft.fftfreq(NX)*NX, ft.fftfreq(NZ)*NZ
#         KX2, KZ2 = KX**2, KZ**2
#         EL, ES = np.zeros(KY.shape), np.zeros(KY.shape)
#         for i in range(NX):
#             for k in range(NZ):
#                 if KX2[i]+KZ2[k]<Kh2:
#                     EL=EL+CS2[i,:,k]
#                 else:
#                     ES=ES+CS2[i,:,k]
                    
#         E=np.sum(CS2,axis=(0,2))
#         spectra[var]['E_large']=EL
#         spectra[var]['E_small']=ES
#         spectra[var]['E_total']=E
        
#         if var=='U1' or var=='U3':
#             spectra[var]['Fr_total']=KY**2*E/Ri_t
#             spectra[var]['Fr_shear']=KY**2*CS2[0,:,0]/Ri_t
        
#         print('Computed spectra for '+var)
#         del CS1,CS2
    
#     return spectra