import numpy as np
import h5py

def c_stats(dir1,dir2):
    fn1, fn2 = dir1+'stats.h5', dir2+'stats.h5'
    f1, f2 = h5py.File(fn1,'r+'), h5py.File(fn2,'r')
    stats = ('/U1rms/','/U2rms/','/U3rms/','/THrms/','/THflux/','/epsilon/','/chi/')
    nk1 = f1[stats[0]].attrs.__getitem__('Samples')
    nk2 = f2[stats[0]].attrs.__getitem__('Samples')
    
    for n in range(1,nk2):
        for stat in stats:
            k=n+nk1-1
            n_old=format(n,"04")
            n_new=format(k,"04")
            f2.copy(stat+n_old, f1[stat], n_new)
            if n==nk2-1:
                f1[stat].attrs.__setitem__('Samples',k+1)
    return()
              
def c_mean(dir1,dir2):
    fn1, fn2 = dir1+'mean.h5', dir2+'mean.h5'
    f1, f2 = h5py.File(fn1,'r+'), h5py.File(fn2,'r')
    means = ('/U1me/','/U3me/','/THme/','/THflux/','/epsilon/','/chi/','/U1U2/','/U3U2/','/U1rms/','/U2rms/','/U3rms/','/THrms/')
    nk1 = f1[means[0]].attrs.__getitem__('Samples')
    nk2 = f2[means[0]].attrs.__getitem__('Samples')
    for n in range(1,nk2):
        for mean in means:
            k=n+nk1-1
            n_old=format(n,"04")
            n_new=format(k,"04")
            f2.copy(mean+n_old, f1[mean], n_new)
            if n==nk2-1:
                f1[mean].attrs.__setitem__('Samples',k+1)
    return()

def c_spectra(dir1,dir2):
    fn1, fn2 = dir1+'spectra.h5', dir2+'spectra.h5'
    f1, f2 = h5py.File(fn1,'r+'), h5py.File(fn2,'r')
    spectra = ('/U1/','/U2/','/U3/','/TH1/')
    nk1 = f1[spectra[0]].attrs.__getitem__('Samples')
    nk2 = f2[spectra[0]].attrs.__getitem__('Samples')
    for n in range(1,nk2):
        for spec in spectra:
            k=n+nk1-1
            n_old=format(n,"04")
            n_new=format(k,"04")
            f2.copy(spec+n_old, f1[spec], n_new)
            if n==nk2-1:
                f1[spec].attrs.__setitem__('Samples',k+1)
    return()

def c_movie(dir1,dir2):
    planes=('xy','xz','yz')
    for plane in planes:
        fn1, fn2 = dir1+'movie_'+plane+'.h5', dir2+'movie_'+plane+'.h5'
        f1, f2 = h5py.File(fn1,'r+'), h5py.File(fn2,'r')
        nk1 = f1.attrs.__getitem__('Samples')
        nk2 = f2.attrs.__getitem__('Samples')
        for n in range(1,nk2):
            k=n+nk1-1
            n_old=format(n,"04")
            n_new=format(k,"04")
            f2.copy('/'+n_old, f1, n_new)
        f1.attrs.__setitem__('Samples',nk1+nk2-1)
    return()