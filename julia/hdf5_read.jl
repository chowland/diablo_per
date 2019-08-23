module h5r

using DelimitedFiles
using HDF5


function read_input(rundir::String)
    fname="$(rundir)input.dat"
    s = readdlm(fname,comments=true, comment_char='C')
    Re=1/s[2,1]
    Ri_t=s[8,1]
    Pr=s[8,2]
    return Re, Ri_t, Pr
end


function read_stat(rundir::String)
    fname="$(rundir)stats.h5"
    h5open(fname,"r") do f
        nk=read(attrs(f["U1rms"]),"Samples")
        tii = zeros(nk)
        U1rms, U2rms, U3rms = zeros(nk), zeros(nk), zeros(nk)
        THrms, THflux = zeros(nk), zeros(nk)
        epsilon, chi = zeros(nk), zeros(nk)
        for i=1:nk
            num=string(i-1,pad=4)
            tii[i]=read(attrs(f["U1rms/$num"]),"Time")
            U1rms[i]=read(f["U1rms/$num"])
            U2rms[i]=read(f["U2rms/$num"])
            U3rms[i]=read(f["U3rms/$num"])
            THrms[i]=read(f["THrms/$num"])
            THflux[i]=read(f["THflux/$num"])
            epsilon[i]=read(f["epsilon/$num"])
            chi[i]=read(f["chi/$num"])
        end
        return tii, U1rms, U2rms, U3rms, THrms, THflux, epsilon, chi
    end
    
end


function read_mean(rundir::String)
    fname="$(rundir)mean.h5"
    h5open(fname,"r") do f
        nk=read(attrs(f["U1me"]),"Samples")
        ny=length(f["U1me/0000"])
        U1me, U3me, THme = zeros(nk,ny), zeros(nk,ny), zeros(nk,ny)
        epsilon, chi, THflux = zeros(nk,ny), zeros(nk,ny), zeros(nk,ny)
        U1U2, U3U2 = zeros(nk,ny), zeros(nk,ny)
        U1U1, U2U2, U3U3, THTH = zeros(nk,ny), zeros(nk,ny), zeros(nk,ny), zeros(nk,ny)
        for i=1:nk
            num=string(i-1,pad=4)
            U1me[i,:]=read(f["U1me/$num"])
            U3me[i,:]=read(f["U3me/$num"])
            THme[i,:]=read(f["THme/$num"])
            epsilon[i,:]=read(f["epsilon/$num"])
            chi[i,:]=read(f["chi/$num"])
            THflux[i,:]=read(f["THflux/$num"])
            U1U2[i,:]=read(f["U1U2/$num"])
            U3U2[i,:]=read(f["U3U2/$num"])
            U1U1[i,:]=read(f["U1rms/$num"])
            U2U2[i,:]=read(f["U2rms/$num"])
            U3U3[i,:]=read(f["U3rms/$num"])
            THTH[i,:]=read(f["THrms/$num"])
        end
        return U1me, U3me, THme, epsilon, chi, THflux, U1U2, U3U2, U1U1, U2U2, U3U3, THTH, nk, ny
    end
    
end




end