% Plot Froude and energy vertical spectra from outxx.h5 files

Ri_t=1;           % Richardson number

fname='../example_run/out01.h5';

dname={'/U1/','/U2/','/U3/','/TH1/'};
S1=h5read(fname,'/U2/');
Kh2=Ri_t/mean(S1(:).^2);  % Horizontal cutoff for large scale spectrum
for n=1:4
    S1=h5read(fname,dname{n});
    [NX,NY,NZ]=size(S1);
    disp(['Loaded variable ' dname{n}])

    % Convert variable to Fourier space
    CS1=fftn(S1/(NX*NY*NZ)); clear S1

    % MATLAB fft indexed (0,1,...,NX/2,-NX/2,...,-2,-1)
    NKY=floor(NY/3);
    KY=0:NKY;

    CS1=0.5*CS1.*conj(CS1); % Calculate 2-sided energy spectrum
    for i=NKY+2:NY-NKY      % Dealias the high wavenumbers
        CS1(i,:,:)=0;
        CS1(:,i,:)=0;
        CS1(:,:,i)=0;
    end
    % Calculate 1-sided energy spectrum:
    CS2=zeros(NX,NKY+1,NZ);

    CS2(:,1,:)=CS1(:,1,:);
    for j=2:NKY+1
        CS2(:,j,:)=CS1(:,j,:)+CS1(:,NY-j+2,:);
    end

    if n==4
        CS2=Ri_t*CS2;
    end
    
    KX=[0:NX/2-1 -NX/2:-1];   KZ=[0:NZ/2-1 -NZ/2:-1];
    E_large{n}=zeros(1,NKY+1); E_small{n}=zeros(1,NKY+1);
    for i=1:NX
        for k=1:NZ
            if KX(i)^2+KZ(k)^2<Kh2
                E_large{n}=E_large{n}+CS2(i,:,k);
            else
                E_small{n}=E_small{n}+CS2(i,:,k);
            end
        end
    end
   
    
    E{n}=sum(sum(CS2,3),1);
    if mod(n,2)==1
        Fr{n}=KY.^2.*E{n}/Ri_t;
        Fr_shear{n}=KY.^2.*CS2(1,:,1)/Ri_t;
    end
    disp(['Computed spectra for ' dname{n}])
    clear CS1 CS2
end

Et_large=E_large{1}+E_large{2}+E_large{3};
Et_small=E_small{1}+E_small{2}+E_small{3};
Et=E{1}+E{2}+E{3};

% figure(1);
% loglog(KY,Fr{1}+Fr{3}); hold on
% loglog(KY,Fr_shear{1}+Fr_shear{3},'--');
% xlabel('$|m|$')
% ylabel('$\Phi_{Fr}$')
% 
% figure(2);
% Et=E{1}+E{2}+E{3};%+E{4};
% loglog(KY,Et.*KY.^3);
% xlabel('$|m|$')
% ylabel('$|m|^3 E(|m|)$')