% Script to combine movie files of consecutive runs
clear
% Target directory where movie file to be extended is
dir1='../run1/';
% Directory containing the movie files for the continued run
dir2='../run2/';

filename ={'movie_xy.h5','movie_xz.h5','movie_yz.h5'}; plane={'z','y','x'};
F={'/U1','/U2','/U3','/TH1'};

f = waitbar(0,'Initializing...','Name','Combining z plane movies');

for idx=1:3     % Loop through the three different plane slices
    
    fname1=[dir1 filename{idx}];
    fname2=[dir2 filename{idx}];
    nk1=h5readatt(fname1,'/','Samples');
    nk2=h5readatt(fname2,'/','Samples');
    
    for n=1:(nk2-1)     % Loop through snapshots from second run
        for i=1:length(F)
            if n<10     % Start reading from 0001 (0000 duplicates end of 1st run)
                dname=['/000' int2str(n) F{i}];
            elseif n<100
                dname=['/00' int2str(n) F{i}];
            elseif n<1000
                dname=['/0' int2str(n) F{i}];
            else
                dname=['/' int2str(n) F{i}];
            end
            G = h5read(fname2,dname);
            [NX,NY]=size(G);
            Time=h5readatt(fname2,dname(1:5),'Time');
            Timestep=h5readatt(fname2,dname(1:5),'Timestep');
            cross=h5readatt(fname2,dname(1:5),plane{idx});
            k=n+nk1-1;
            if i==1
                waitbar(double(n)/double(nk2-1)/3+double(idx-1)/3,f,['Writing time step ' int2str(k)])
            end
            if k<10
                dname=['/000' int2str(k) F{i}];
            elseif k<100
                dname=['/00' int2str(k) F{i}];
            elseif k<1000
                dname=['/0' int2str(k) F{i}];
            else
                dname=['/' int2str(k) F{i}];
            end
            h5create(fname1,dname,[NX NY],'ChunkSize',[NX 1]);
            h5write(fname1,dname,G);
            h5writeatt(fname1,dname(1:5),'Time',Time);
            h5writeatt(fname1,dname(1:5),'Timestep',Timestep);
            h5writeatt(fname1,dname(1:5),plane{idx},cross);
        end
    end
    h5writeatt(fname1,'/','Samples',nk1+nk2-1);
    if idx<3
        waitbar(double(idx)/3,f,'Name',['Combining ' plane{idx+1} ' plane movies']);
    end
end
delete(f);