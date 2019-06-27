import numpy as np
import h5py

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

def read_stat(rundir):
    fname=rundir+'stats.h5'
    f=h5py.File(fname,'r')
    nk=f['/U1rms'].attrs.__getitem__('Samples')
    tii=np.zeros(nk)
    U1rms, U2rms, U3rms = np.zeros(nk), np.zeros(nk), np.zeros(nk)
    THrms, THflux = np.zeros(nk), np.zeros(nk)
    epsilon, chi = np.zeros(nk), np.zeros(nk)
    for i in range(nk):
        if i<10:
            num='000'+str(i)
        elif i<100:
            num='00'+str(i)
        elif i<1000:
            num='0'+str(i)
        else:
            num=str(i)
        tii[i]=f['/U1rms/'+num].attrs.__getitem__('Time')
        U1rms[i]=f['/U1rms/'+num][()]
        U2rms[i]=f['/U2rms/'+num][()]
        U3rms[i]=f['/U3rms/'+num][()]
        THrms[i]=f['/THrms/'+num][()]
        THflux[i]=f['/THflux/'+num][()]
        epsilon[i]=f['/epsilon/'+num][()]
        chi[i]=f['/chi/'+num][()]
    return (tii, U1rms, U2rms, U3rms, THrms, THflux, epsilon, chi)

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
        if i<10:
            num='000'+str(i)
        elif i<100:
            num='00'+str(i)
        elif i<1000:
            num='0'+str(i)
        else:
            num=str(i)
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

def read_spectra(rundir):
    fname=rundir+'spectra.h5'
    f=h5py.File(fname,'r')
    nk=f['/U1'].attrs.__getitem__('Samples')
    nky=(f['/U1/0000'].size-1)/2
    U1, U2 = np.zeros((nk,2*nky+1)), np.zeros((nk,2*nky+1))
    U3, TH = np.zeros((nk,2*nky+1)), np.zeros((nk,2*nky+1))
    for i in range(nk):
        if i<10:
            num='000'+str(i)
        elif i<100:
            num='00'+str(i)
        elif i<1000:
            num='0'+str(i)
        else:
            num=str(i)
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
