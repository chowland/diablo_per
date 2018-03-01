% Read movie slices from movie_ij.h5

% Run directory
rundir = '/local/scratch/public/cjh225/CRe1e4_1e3/';

% Movie file (choose plane)
plane='xy';
% Quantity to plot ('/U', '/V', '/W' or '/THn')
F = '/W';

fname = [rundir 'movie_' plane '.h5'];
% Number of samples
nk=h5readatt(fname,'/','Samples');
if strcmp(plane,'xy')
    c_plane = h5readatt(fname,'/0000','z');
elseif strcmp(plane,'xz')
    c_plane = h5readatt(fname,'/0000','y');
elseif strcmp(plane,'yz')
    c_plane = h5readatt(fname,'/0000','x');
end

for i=1:nk
    if i<=10
        dname=['/000' int2str(i-1) F];
    elseif i<=100
        dname=['/00' int2str(i-1) F];
    elseif i<=1000
        dname=['/0' int2str(i-1) F];
    else
        dname=['/' int2str(i-1) F];
    end
    G = h5read(fname,dname);
    xvec = linspace(0,2*pi,size(G,1));
    yvec = linspace(0,2*pi,size(G,2));
    if strcmp(plane,'yz')
        G = G(:,:)';
    end
    pcolor(xvec,yvec,G(:,:)'); shading interp
    if i==1
        c = max(abs(G(:)));
    end
    caxis([-c c])
    if strcmp(F,'/TH1');
        colormap(cmocean('-dense'))
        title('$\theta$');
    else
        colormap(cmocean('balance'))
        if strcmp(F,'/U')
            title('$u$')
        elseif strcmp(F,'/V')
            title('$v$')
        elseif strcmp(F,'/W')
            title('$w$')
        end
    end
    colorbar
    if strcmp(plane,'xy')
        xlabel('$x$'); ylabel('$y$');
    elseif strcmp(plane,'xz')
        xlabel('$x$'); ylabel('$z$');
    elseif strcmp(plane,'yz')
        xlabel('$z$'); ylabel('$y$');
    end
    M(i) = getframe;
    if i~=nk
        clf;
    end
end