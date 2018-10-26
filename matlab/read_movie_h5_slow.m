% Read movie slices from movie_ij.h5

rundir='../example_run/';

% Select plane and variable to read (U1, U2, U3, TH1)
plane='xy';     F='/U1';

fname=[rundir 'movie_' plane '.h5'];

nk=h5readatt(fname,'/','Samples');

for n=1:nk
    % Create the right dataset name for each time step
    if n<=10
        dname=['/000' int2str(n-1) F];
    elseif n<=100
        dname=['/00' int2str(n-1) F];
    elseif n<=1000
        dname=['/0' int2str(n-1) F];
    else
        dname=['/' int2str(n-1) F];
    end

    % Read slices & transpose if necessary
    G=h5read(fname,dname);
    if strcmp(plane,'yz');  G = G'; end
    [NX,NY]=size(G);
    if n==1
        yvec = linspace(2*pi/NY,2*pi,NY);
    end

    % Add background stratification to scalar field (optional)
    if strcmp(F,'/TH1') && ~strcmp(plane,'xz')
        for j=1:length(yvec)
            G(:,j)=G(:,j)+yvec(j);
        end
    end

    pcolor(yvec,yvec,G'); shading flat; colorbar
    
    if strcmp(F,'/TH1')
        if strcmp(plane,'xz')
            caxis([-0.5 0.5])
        else
            caxis([0 2*pi])
        end
        cmocean('-dense');
    else
        caxis([-1 1])
        cmocean('balance');
    end
    M(n)=getframe;
end