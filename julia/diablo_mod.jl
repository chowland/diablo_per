module out_analysis

using HDF5
using FFTW
using Statistics


function out2chi(fname)
    
    # Read variable from out*.h5 file
    f = h5open(fname,"r")
    (NX,NY,NZ) = size(f["TH1"])
    S1=read(f["TH1"])
    close(f)
    
    # Subtract horizontal mean from variable
    Sme=mean(S1,dims=(1,3))
    S1=S1.-Sme
    
    # Preallocate chi variable
    chi=zeros(NX,NY,NZ)

    # Create FFT plans
    FFTW.set_num_threads(32)
    ft1=plan_rfft(S1,1)
    ft2=plan_rfft(S1,2)
    ft3=plan_rfft(S1,3)
    
    # Compute derivatives and add to chi
    CS1=ft1 * S1
    CIKX=im*(0:NX÷2)
    CS1=CIKX.*CS1
    chi+=(inv(ft1) * CS1).^2
    CS1=nothing
    println("Computed 1st derivative")

    CS1=ft2 * S1
    CIKY=reshape(im*(0:NY÷2),(1,NY÷2+1,1))
    CS1=CIKY.*CS1
    chi+=(inv(ft2) * CS1).^2
    CS1=nothing
    println("Computed 2nd derivative")

    CS1=ft3 * S1
    CIKZ=reshape(im*(0:NZ÷2),(1,1,NZ÷2+1))
    CS1=CIKZ.*CS1
    chi+=(inv(ft3) * CS1).^2
    CS1=nothing
    println("Computed 3rd derivative")

    Re, Ri_t, Pr = 1e4, 1, 1
    chi=Ri_t/Re/Pr*chi

    # Clear large temporary variable
    CS1=nothing
    
    return chi
end

function out2eps(fname)
    
    dname=("U1","U2","U3")
    
    for i=1:3
        # Read variable from out*.h5 file
        f = h5open(fname,"r")
        if i==1
            (NX,NY,NZ) = size(f[dname[i]])
        end
        S1=read(f[dname[i]])
        close(f)
    
        # Subtract horizontal mean from variable
        if mod(i,2)==1
            Sme=mean(S1,dims=(1,3))
            S1=S1.-Sme
        end
    
        # Preallocate chi variable
        if i==1
            epsilon=zeros(NX,NY,NZ)
        end

        # Create FFT plans
        if i==1
            FFTW.set_num_threads(32)
            ft1=plan_rfft(S1,1)
            ft2=plan_rfft(S1,2)
            ft3=plan_rfft(S1,3)
        end
    
        # Compute derivatives and add to chi
        CS1=ft1 * S1
        println("Performed ft1")
        CIKX=im*(0:NX÷2)
        println("Created ik vector")
        CS1=CIKX.*CS1
        println("Multiplied ik")
        epsilon=epsilon+(inv(ft1) * CS1).^2
        println("Computed 1st derivative")
        CS1=nothing

        CS1=ft2 * S1
        println("Performed ft2")
        CIKY=reshape(im*(0:NY÷2),(1,NY÷2+1,1))
        println("Created ik vector")
        CS1=CIKY.*CS1
        println("Multiplied ik")
        epsilon=epsilon+(inv(ft2) * CS1).^2
        CS1=nothing

        println("Computed 2nd derivative")

        CS1=ft3 * S1
        CIKZ=reshape(im*(0:NZ÷2),(1,1,NZ÷2+1))
        CS1=CIKZ.*CS1
        epsilon=epsilon+(inv(ft3) * CS1).^2
        CS1=nothing

        println("Computed 3rd derivative")
    
        # Clear large temporary variable
        CS1=nothing
    
    end

    Re = 1e4
    epsilon=1/Re*epsilon

    return epsilon
end

end