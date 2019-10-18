import numpy as np
import matplotlib.image as img
import h5py
import cmocean

def make_movie(rundir,plane,var):
    fname=rundir+'movie_'+plane+'.h5'
    f=h5py.File(fname,'r')
    nk=f.attrs.__getitem__('Samples')
    for i in range(nk):
        dname='/'+format(i,"04")+"/"+var
        G=f[dname][()]
        
        if plane=='yz':
            G=G.T
        if i==0:
            NX,NY=G.shape
            yvec=np.linspace(2*np.pi/NY,2*np.pi,NY)
            
        if var=='TH1' and plane!='xz':
            for j in range(NY):
                G[j,:]=G[j,:]+yvec[j]
        
        if var=='TH1':
            if plane=='xz':
                if i==0:
                    c=0.4
                im=np.flipud(G)
                img.imsave('tmp/'+dname[1:5]+'.png',im,vmin=-c,vmax=c,cmap=cmocean.cm.deep_r)
            else:
                if i==0:
                    c=2*np.pi
                #img.imsave('tmp/'+dname[1:5]+'.png',G,vmin=0,vmax=c,cmap=cmocean.cm.ice,origin='lower')
                img.imsave('tmp/'+dname[1:5]+'.png',G,vmin=np.min(G),vmax=np.max(G),cmap=cmocean.cm.deep_r,origin='lower')
        elif var=='U2':
            if i==0:
                c=0.4
            img.imsave('tmp/'+dname[1:5]+'.png',G,vmin=-c,vmax=c,cmap=cmocean.cm.balance,origin='lower')
        else:
            if i==0:
                c=1
            img.imsave('tmp/'+dname[1:5]+'.png',G,vmin=-c,vmax=c,cmap=cmocean.cm.balance,origin='lower')     
            
import pyfftw
            
def vorticity_movie(rundir):
    fname=rundir+'movie_xy.h5'
    f=h5py.File(fname,'r')
    nk=f.attrs.__getitem__('Samples')
    NX, NY = 1024, 1024
    
    A1=pyfftw.empty_aligned((NX,NY),dtype='float64')
    A2=pyfftw.empty_aligned((NX,NY),dtype='float64')
    B1=pyfftw.empty_aligned((NX//2+1,NY),dtype='complex128')
    B2=pyfftw.empty_aligned((NX,NY//2+1),dtype='complex128')
    
    ft1=pyfftw.FFTW(A1,B1,axes=(0,), direction='FFTW_FORWARD',threads=4)
    ift1=pyfftw.FFTW(B1,A1,axes=(0,), direction='FFTW_BACKWARD',threads=4)
    ft2=pyfftw.FFTW(A2,B2,axes=(1,), direction='FFTW_FORWARD',threads=4)
    ift2=pyfftw.FFTW(B2,A2,axes=(1,), direction='FFTW_BACKWARD',threads=4)
    
    CIKX=1j*np.arange(NX//2+1).reshape(NX//2+1,1)
    CIKY=1j*np.arange(NY//2+1).reshape(1,NY//2+1)
    
    for i in range(nk):
        dname='/'+format(i,"04")+"/U1"
        U1=f[dname][()].T
        dname='/'+format(i,"04")+"/U2"
        U2=f[dname][()].T
        
        CS1=ft1(U2)
        CS1=CIKX*CS1
        Vx=ift1(CS1)
        
        CS1=ft2(U1)
        CS1=CIKY*CS1
        Uy=ift2(CS1)
        
        Omega=Vx-Uy
        
        img.imsave('tmp/'+dname[1:5]+'.png',Omega.T,
                   vmin=-3,vmax=3,cmap=cmocean.cm.curl,origin='lower')
        if i%10==0:
            print('Saved frame '+str(i))
            
def thorpe_2d(rundir):
    fname=rundir+'movie_xy.h5'
    f=h5py.File(fname,'r')
    nk=f.attrs.__getitem__('Samples')
    L_T2=np.zeros(nk)
    for i in range(nk):
        dname='/'+format(i,"04")+'/TH1/'
        TH=f[dname][()]
        
        if i==0:
            NX,NY=TH.shape
            yvec=np.linspace(2*np.pi/NY,2*np.pi,NY)
            
        for j in range(NY):
            TH[j,:]=TH[j,:]+yvec[j]
            
        delta=np.zeros((NX,NY))
        THsort=TH.flatten()
        THsort.sort()
        for k in range(NX):
            for j in range(NY):
                n=np.mean(np.nonzero(THsort==TH[k,j]))
                delta[k,j]=abs(n/np.float(NX)-j)
                
        delta=delta*2*np.pi/NY
        L_T2[i]=np.sqrt(np.mean(delta**2))
        print('Time '+str(i)+' complete')
    return L_T2

def thorpe_1d(rundir,idx):
    fname=rundir+'movie_xy.h5'
    f=h5py.File(fname,'r')
    nk=f.attrs.__getitem__('Samples')
    L_T=np.zeros(nk)
    for i in range(nk):
        dname='/'+format(i,"04")+'/TH1/'
        TH=f[dname][()]
        
        if i==0:
            NX,NY=TH.shape
            yvec=np.linspace(2*np.pi/NY,2*np.pi,NY)
            
        for j in range(NY):
            TH[j,:]=TH[j,:]+yvec[j]
            
        THsort=np.argsort(TH[:,idx])
        delta=THsort-np.arange(NY)
                
        delta=delta*2*np.pi/NY
        L_T[i]=np.sqrt(np.mean(delta**2))
    return L_T

def thorpe_profile(rundir,idx,n_time):
    fname=rundir+'movie_xy.h5'
    f=h5py.File(fname,'r')
    nk=f.attrs.__getitem__('Samples')
    L_T=np.zeros(nk)
    dname='/'+format(i,"04")+'/TH1/'
    TH=f[dname][()]

    NX,NY=TH.shape
    yvec=np.linspace(2*np.pi/NY,2*np.pi,NY)
            
    for j in range(NY):
        TH[j,:]=TH[j,:]+yvec[j]
            
    THsort=np.argsort(TH[:,idx])
    THog=TH[:,idx]
    delta=THsort-np.arange(NY)
    THstd=np.zeros(np.size(THog))
    for j in range(NY):
        THstd[j]=THog[THsort[j]]
                
    L_T=delta*2*np.pi/NY
    return L_T, THstd, THog