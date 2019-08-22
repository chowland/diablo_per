% Plot horizontal spectra

fname='../example_run/end.h5';

dname={'/U1/','/U2/','/U3/','/TH1/'};

for n=1:4
    S1=h5read(fname,dname{n});
    [NX,NY,NZ]=size(S1);
    disp(['Loaded variable ' dname{n}])
    
    % Convert variable to Fourier space
    CS1=fftn(S1/(NX*NY*NZ)); clear S1

    % MATLAB fft indexed (0,1,...,NX/2,-NX/2,...,-2,-1)
    NK=floor(NX/3);
    KH=0:NK;

    CS1=0.5*CS1.*conj(CS1);
    CS1(NK+2:NX-NK,:,:)=0;
    CS1(:,NK+2:NY-NK,:)=0;
    CS1(:,:,NK+2:NZ-NK)=0;

    KX=[0:NX/2-1 -NX/2:-1];     KZ=[0:NZ/2-1 -NZ/2:-1];
    KX2=KX.^2;                  KZ2=KZ.^2;

    CS2=zeros(NK+1,NY);
    for i=1:NX
        for k=1:NZ
            m=round(sqrt(KX2(i)+KZ2(k)));
            if m<=NK
                CS2(m+1,:)=CS2(m+1,:)+CS1(i,:,k);
            end
        end
    end

    Eh{n}=sum(CS2,2);
    disp(['Computed spectra for ' dname{n}])
    clear CS1 CS2
end

Eht=Eh{1}+Eh{2}+Eh{3};
figure;
loglog(KH,Eht,'k'); hold on
loglog(KH,Eh{1},'r--')
loglog(KH,Eh{2},'g--')
loglog(KH,Eh{3},'b--')