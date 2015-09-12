function [cgSoln, cfSoln, cgOn, cfOn] = PDELoadSolnFiles(fDir, onflag)
% PDELoadSolnFiles: Loads in data for multiple instances of the optimal stopping
% problem. Calls PDESolnLoad() for each instance to be loaded. Typically
% called after PDECreateSolnFiles() is called (which creates all the file
% to be loaded).
% INPUTS: 
%   fDir: base name for directory to hold the files with the solutions to the several PDEs.
%   onflag: if true, then online learning solutions are also created, if
%       false, then only offline learning solutions are created.
%
% OUTPUTS:
%   CGOff: PDE solution structure for offline learning, 0 sampling cost,
%      positive discount rate
%   CFOff: PDE solution structure for offline learning, positive sampling cost,
%      zero discount rate
%   CGOn: PDE solution structure for ONline learning, 0 sampling cost,
%      positive discount rate (or [] if onflag is false)
%   CFOn: PDE solution structure for ONline learning, positive sampling cost,
%      zero discount rate (or [] if onflag is false)
%   Side effect: a bunch of plots as calculations are run - if they don't
%       'look right' then there is probably a numerical stability issue to
%       address.
%
% Code provided 'as is' with no warrantees of correctness.
%
% Author: S Chick, 2015, all rights reserved.
% Created: 2015 08 17.

PDELocalInit;

if nargin < 1
    fDir = PDEmatfilebase;
end
if nargin < 2
    onflag = false;  % default is to only generate files for offline learning
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First do offline files
[rval, cfSoln] = PDESolnLoad([fDir PDEnodiscbase PDEoffbase]); % Load in the undiscounted reward offline case
[rval, cgSoln] = PDESolnLoad([fDir PDEdiscbase PDEoffbase]); % Load in the discounted reward offline case

% Next do online files
if onflag
    [rval, cfOn] = PDESolnLoad([fDir PDEnodiscbase PDEonbase]); % Load in the undiscounted reward online case
    [rval, cgOn] = PDESolnLoad([fDir PDEdiscbase PDEonbase]); % Load in the discounted reward online case
else
    cfOn = [];
    cgOn = [];
end

end
