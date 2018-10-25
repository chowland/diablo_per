% Read parameters in from input files

fname=[rundir 'input.dat'];
fileID=fopen(fname,'r');
formatSpec='%f%f%f%f';
A=textscan(fileID,formatSpec,1,'Delimiter',' ', ...
    'MultipleDelimsAsOne',true,'HeaderLines',6);
fclose(fileID);
NU=A{1}; clear A

fileID=fopen(fname,'r');
formatSpec='%f%f%f';
A=textscan(fileID,formatSpec,1,'Delimiter',' ', ...
    'MultipleDelimsAsOne',true,'HeaderLines',19);
fclose(fileID);
Ri_t=A{1};     Pr=A{2}; clear A

clear fname fileId formatSpec