function [ rval ] = PDELocalSetPaths( basepath )
% PDELocalSetPaths: 
%   PURPOSE: updates the matlab path variable so that the code for the
%   solution of free boundary problems of the sort delt with to approximate
%   the optimal solution to certain sequential sampling problems.
%
%   INPUT: basepath is either empty (to default to the current directory in
%   the matlab environment) or a text string with the path name associated
%   with the routines. Do not terminate the directory name with a \ on
%   windows systems.
%   OUTPUT: always returns true
%   EXAMPLE USAGE:
%       LocalDelaySetPaths();
%       LocalDelaySetPaths('d:\users\localforkofrepo');
%   INTENDED WORK FLOW:
%       When downloading matlab code for the paper, this file should be
%       updated once to adapt to the local machine's requirements. This
%       update is not to be committed back to the code repo.
%
% Source provided 'as is' with no warrantees or claims provided or implied.
% 2015 S Chick

    if nargin < 1
        BASEDIR = pwd;  
    else
        BASEDIR = basepath; % 
    end

    delaycoredir = [ BASEDIR '\pdecode\'];
    addpath(genpath(delaycoredir));

    delaypaperdir = [ BASEDIR '\papercg\'];
    addpath(genpath(delaypaperdir));

    delaypaperdir = [ BASEDIR '\papercf\'];
    addpath(genpath(delaypaperdir));

    delaypaperdir = [ BASEDIR '\papercfp\'];
    addpath(genpath(delaypaperdir));

    basedelaydir = [BASEDIR];    
    addpath(basedelaydir);              % set up this way to preclude .git directory structure from being added to path
    
    rval = 1;
end