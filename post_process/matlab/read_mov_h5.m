% Run directory
rundir = '/local/scratch/public/cjh225/plane_wave/';
% Movie file (choose plane)
fname = [rundir 'movie_xz.h5'];
% Quantity to plot ('/U', '/V', '/W' or '/THn')
F = '/W';

% Number of samples
nk=h5readatt(fname,'/','Samples');
plane = h5readatt(fname,'/0000','y');

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
    pcolor(xvec,yvec,G(:,:)'); shading interp
%     if i==1
        c=max(abs(G(:)));
%     end
    caxis([-c c])
    colormap(cmocean('balance'))
    colorbar
%     title(['$z = ' num2str(plane) ',\, t = ' num2str(tii(i),'%05.1f') '$'])
    M(i) = getframe;
    if i~=nk
        clf;
    end
end