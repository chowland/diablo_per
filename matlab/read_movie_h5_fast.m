% Read movie slices from movie_ij.h5

rundir='../example_run/';

% Select plane and variable to read (U1, U2, U3, TH1) & whether to write to
% AVI or not (edit parameters for video at end of script)
plane='xy';     F='/U1';    write_out=true;

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

    if strcmp(F,'/TH1')
        if strcmp(plane,'xz')
            c=1/2;      % Set colour axis to go from -c to c
            im=flipud((G'+c)*256/(2*c)+1);
        else
            c = 2*pi;   % Set colour axis to go from 0 to c
            im=flipud(G'*256/c);
        end
        for i=1:NX
            for j=1:NY
                if im(i,j)<1
                    im(i,j)=1;
                elseif im(i,j)>256
                    im(i,j)=256;
                end
            end
        end
        mov(:,:,n)=im;
        M(n)=im2frame(im,cmocean('-dense'));

    else
        if n==1;    c=1;    end     % Set colour axis to go from -c to c
        im=flipud((G'+c)*256/(2*c)+1);
        for i=1:NX
            for j=1:NY
                if im(i,j)<1
                    im(i,j)=1;
                elseif im(i,j)>256
                    im(i,j)=256;
                end
            end
        end
        mov(:,:,n)=im;
        M(n)=im2frame(im,cmocean('balance'));
    end
    
end

% Optionally write to an AVI file
if write_out
    filename=[F(2:end) '_' plane '.avi'];
    myVideo=VideoWriter(filename);
    myVideo.Quality=100;
    myVideo.FrameRate=15;
    open(myVideo)
    writeVideo(myVideo,M)
    close(myVideo)
end