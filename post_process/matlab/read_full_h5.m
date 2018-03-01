% Read outxx.h5 files containing full flow fields

% Run directory
rundir = '/local/scratch/public/cjh225/BRe1e4_1e3/';
% rundir = '../../../scratch_backup/large_force/';

fname = [rundir 'restart_files/out50.h5'];
% fname = [rundir 'start.h5'];

U = h5read(fname,'/U/');
V = h5read(fname,'/V/');
W = h5read(fname,'/W/');
TH= h5read(fname,'/TH1/');