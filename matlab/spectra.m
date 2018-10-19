% Plot Froude and energy vertical spectra from outxx.h5 files

Ri_t=1;           % Richardson number

fname='../example_run/out01.h5';
dname={'/U1/','/U2/','/U3/','/TH1/'};
for n=1:4
    S1=h5read(fname,dname{n});
    [NX,NY,NZ]=size(S1);
    disp(['Loaded variable ' dname{n}])

    % Convert variable to Fourier space
    CS1=fftn(S1/(NX*NY*NZ)); clear S1

    % MATLAB fft indexed (0,1,...,NX/2,-NX/2,...,-2,-1)
    NKY=floor(NZ/3);
    KY=0:NKY;
    Kh=35;  % Horizontal cutoff for large scale spectrum

    CS1=0.5*CS1.*conj(CS1); % Calculate 2-sided energy spectrum
    
    % Calculate 1-sided energy spectrum:
    CS2=zeros(NX,NKY+1,NZ);

    CS2(:,1,:)=CS1(:,1,:)
    for j=2:NKY+1
        CS2(:,j,:)=CS1(:,j,:)+CS1(:,NY-j+2,:);
    end

    if n==4
        CS2=RI*CS2;
    end
    E{n}=sum(sum(CS2,3),1);
    if n<=2
        Fr{n}=KY.^2.*E{n}/RI;
        Fr_shear{n}=KY.^2.*CS2(1,:,1)/RI;
    end
    disp(['Computed spectra for ' dname{n}])
    clear CS1 CS2
end

figure(1);
loglog(KY,Fr{1}+Fr{2}); hold on
loglog(KY,Fr_shear{1}+Fr_shear{2},'--');
xlabel('$|m|$')
ylabel('$\Phi_{Fr}$')

figure(2);
Et=E{1}+E{2}+E{3}+E{4};
loglog(KY,Et.*KY.^3);
xlabel('$|m|$')
ylabel('$|m|^3 E(|m|)$')