
% Run as a macro in order to set up variables which are used in a variety
% of routines.
PDEdiscbase = 'CG';
PDEnodiscbase = 'CF';
PDEonbase = 'On';
PDEoffbase = '';
%[basefolder,name,ext] = fileparts(pwd);
[basefoldermfile,name,ext] = fileparts(mfilename('fullpath'));
basefolder = fileparts(basefoldermfile);
PDEmatfilebase = [ basefolder '\Matfiles\'];

PDEfigfilebase = 'Figure\';

addpath(genpath(basefolder));
