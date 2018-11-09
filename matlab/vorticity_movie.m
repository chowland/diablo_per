% Plot movie of vorticity component of a plane slice

% Run directory
rundir='../example_run/';

% Select plane to read (xy, xz or yz) to read & whether to write to AVI or
% not (edit paratemers for video at end of script)
plane='xz';     write_out=true;

fname = [rundir 'movie_' plane '.h5'];

nk=h5readatt(fname,'/','Samples');

if strcmp(plane,'xy')
    c_plane = h5readatt(fname,'/0000','z');
    F1='/U2'; F2='/U1';
elseif strcmp(plane,'xz')
    c_plane = h5readatt(fname,'/0000','y');
    F1='/U1'; F2='/U3';
elseif strcmp(plane,'yz')
    c_plane = h5readatt(fname,'/0000','x');
    F1='/U3'; F2='/U2';
end

f = waitbar(0,'Initializing...','Name','Creating vorticity movie');

for i=1:nk
    if i<=10
        dname1=['/000' int2str(i-1) F1];
        dname2=['/000' int2str(i-1) F2];
    elseif i<=100
        dname1=['/00' int2str(i-1) F1];
        dname2=['/00' int2str(i-1) F2];
    elseif i<=1000
        dname1=['/0' int2str(i-1) F1];
        dname2=['/0' int2str(i-1) F2];
    else
        dname1=['/' int2str(i-1) F1];
        dname2=['/' int2str(i-1) F2];
    end
    G1 = h5read(fname,dname1);
    G2 = h5read(fname,dname2);
    if strcmp(plane,'yz')
        G1 = G1';  G2 = G2';
    end
    if i==1
        [NX,NY]=size(G1);
        xvec = linspace(0,2*pi,NX);
        yvec = linspace(0,2*pi,NY);
        CIKX = 1i*[0:ceil((NX-1)/2) ceil(-(NX-1)/2):-1];
        CIKY = 1i*[0:ceil((NY-1)/2) ceil(-(NY-1)/2):-1];
    end
    CG1=fftn(G1);   CG2=fftn(G2);
    for j=1:NX
        for k=1:NY
            C_Omega(j,k)=CG1(j,k)*CIKY(k)-CG2(j,k)*CIKX(j);
        end
    end
    Omega=ifftn(C_Omega,'symmetric');
    if strcmp(plane,'xy')
        Omega=-Omega;
    end

    c = 10;
    
    im=flipud((Omega'+c)*256/(2*c)+1);
    for k=1:NX
        for j=1:NY
            if im(k,j)<1
                im(k,j)=1;
            elseif im(k,j)>256
                im(k,j)=256;
            end
        end
    end
    M(i)=im2frame(im,cmocean('curl'));
    
    waitbar(double(i)/double(nk),f,['Calculating frame ' int2str(i) ' of ' int2str(nk)])
    if i==nk
        delete(f)
    end

end

if write_out
    if strcmp(plane,'yz')
        filename='omega_x.avi';
    elseif strcmp(plane,'xy')
        filename='omega_y.avi';
    else
        filename='omega_z.avi';
    end
    myVideo=VideoWriter(filename);
    myVideo.Quality=100;
    myVideo.FrameRate=15;
    open(myVideo)
    writeVideo(myVideo,M)
    close(myVideo)
end