import numpy as np
import matplotlib.image as img
import h5py
import cmocean

def make_movie(rundir,plane,var):
    fname=rundir+'movie_'+plane+'.h5'
    f=h5py.File(fname,'r')
    nk=f.attrs.__getitem__('Samples')
    for i in range(nk):
        if i<10:
            dname='/000'+str(i)+'/'+var
        elif i<100:
            dname='/00'+str(i)+'/'+var
        elif i<1000:
            dname='/0'+str(i)+'/'+var
        else:
            dname='/'+str(i)+'/'+var
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
                    c=0.5
                im=np.flipud(G)
                img.imsave('tmp/'+dname[1:5]+'.png',im,vmin=-c,vmax=c,cmap=cmocean.cm.ice)
            else:
                if i==0:
                    c=2*np.pi
                img.imsave('tmp/'+dname[1:5]+'.png',G,vmin=0,vmax=c,cmap=cmocean.cm.ice,origin='lower')
        else:
            if i==0:
                c=1
            img.imsave('tmp/'+dname[1:5]+'.png',G,vmin=-c,vmax=c,cmap=cmocean.cm.balance,origin='lower')
