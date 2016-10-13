function [rval] = PDECreateSolnFiles( fDir, onflag, THoriz)
% PDECreateSolnFiles: Solves free boundary problem in scaled coordinates for the
% sequential sampling for simulation optimization. 
%
% For sample calling conventions, including how to set up the fields of the
% structures PDEscale and PDEparam, please see TestSolvePDE.m.
%
% INPUTS: 
%   fDir: base name for directory to hold the files with the solutions to the several PDEs.
%   onflag: if true, then online learning solutions are also created, if
%       false, then only offline learning solutions are created.
%   THoriz: set to positive, finite value if finite time horizon is desired, 
%       with given number of replications. If 0 or negative, then infinite 
%       time horizon is assumed. THoriz set assuming params scaled so that
%       sigma = 1, offline learning, and (c=0, disc=10^-3) for discounted
%       and (c=1,disc=0) for undiscounted.
%
% OUTPUTS:
%   Side effect: files of name fName[On/Off]<n>.mat for n=1,2,...,numfiles
%       which contain solutions for PDE in standardized coordinates 
%       Contents include: upper and lower boundaries, standardized: reward
%       function, PCS, expected number of samples to stop time, ...
%       also creates fName[On/Off]0.mat, with summary information for files.
%           base name (with On/Off) can be found in <structurename>.Header.fName
%
% Online learning: not yet supported.
% CFP: not yet supported.
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

if nargin < 3
    THoriz = -1;    % default is to run infinite horizon analysis
end

if THoriz <= 0;
    finiteT = false; % by default, compute infinite horizon solution
else
    finiteT = true;
end

PDELocalInit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up generic functions for zero discounting,
% these functions should return two vectors of the same size as wvec: the
% first should be the expected reward at time s given wvec, the second
% should be the expected number of samples to take additionally to achieve
% it.
% First set things up for offline learning, then do online learning files
% if needed.

generictermreward=@(wvec,s,p1,p2) PDEsimplereward(wvec);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
baseparams = { 'online', 0, 'retire', 0, 'DoPlot', 1 };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                               %% 
%%       Create functions for use with zero discounting          %%   offline case
%%                                                               %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CFApproxValuefunc=@(wvec,s,p1,p2) PDECFApproxValue(wvec,s,p1);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
upperNoDisc=@(s,p1,p2) CFApproxBoundW(s);
CFscalevec = {'c', 1, 'sigma', 10e5, 'discrate', 0, 'P', 1};
if finiteT
    CFfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', generictermreward, 'approxmethod', upperNoDisc}; % use this to not use KG* for terminal reward at time 'infinity'
    CFparamvec = { 't0', .1, 'tEND', THoriz, 'precfactor', 10, 'ceilfactor', 1.1, 'BaseFileName', [PDEnodiscbase PDEoffbase] , 'matdir' , PDEmatfilebase, 'finiteT', finiteT  };
else
    CFfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', CFApproxValuefunc, 'approxmethod', upperNoDisc}; % use this to have KG* type rule at time 'infinity' for ca
%    CFparamvec = { 't0', .1, 'tEND', 100000, 'precfactor', 10, 'ceilfactor', 1.1, 'BaseFileName', [PDEnodiscbase PDEoffbase] , 'matdir' , PDEmatfilebase, 'finiteT', finiteT  };
    CFparamvec = { 't0', .1, 'tEND', 100000, 'precfactor', 12, 'ceilfactor', 1.15, 'BaseFileName', [PDEnodiscbase PDEoffbase] , 'matdir' , PDEmatfilebase, 'finiteT', finiteT  };
    %figdir Figure\,  UnkVariance 0
end
scalevec = CFscalevec; 
paramvec = [CFparamvec, CFfunctionset, baseparams];
[scale, param] = PDEInputConstructor( scalevec, paramvec );
[scale, param] = PDEInputValidator( scale, param )
s0 = 1/(scale.gamma * param.t0)
sEND = 1/(scale.gamma * param.tEND)
tic
[rval, MAXFiles] = PDESolnCompute(scale, param);
% Load in the data structures form those computations
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                               %% 
%%     Create functions for use with positive discounting        %%  offline case
%%                                                               %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up generic functions for positive discounting
CGApproxValuefunc=@(wvec,s,p1,p2) PDECGApproxValue(wvec,s,p1);   % this is valid terminal reward for discounted rewards, valued in time s currency
upperDisc=@(s,p1,p2) CGApproxBoundW(s);
CGscalevec = {'c', 0, 'sigma', 10e5, 'discrate', 0.0002, 'P', 1 };
if finiteT
    CGfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', generictermreward, 'approxmethod', upperDisc};
    CGparamvec = { 't0', 0.002, 'tEND', THoriz, 'precfactor', 8, 'ceilfactor', 1.8, 'BaseFileName', [PDEdiscbase PDEoffbase] , 'matdir' , PDEmatfilebase};
else
    CGfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', CGApproxValuefunc, 'approxmethod', upperDisc};
    CGparamvec = { 't0', 0.002, 'tEND', 400000, 'precfactor', 8, 'ceilfactor', 1.8, 'BaseFileName', [PDEdiscbase PDEoffbase] , 'matdir' , PDEmatfilebase};
    %figdir Figure\,  UnkVariance 0
end

scalevec = CGscalevec; 
paramvec = [CGparamvec, CGfunctionset, baseparams];
[scale, param] = PDEInputConstructor( scalevec, paramvec );
[scale, param] = PDEInputValidator( scale, param )
tic
[rval, MAXFiles] = PDESolnCompute(scale, param);
toc

%%%%%%%%%%%%%%

if onflag   % online files are requested
	% do first online with zero discounting
    % Online with zero discounting is not ok unless time horizon is finite
    CFApproxValuefunc=@(wvec,s,p1,p2) PDECFApproxValue(wvec,s,p1);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
    upperNoDisc=@(s,p1,p2) CFApproxBoundW(s);
    CFscalevec = {'c', 1, 'sigma', 10e5, 'discrate', 0, 'P', 1};
	if finiteT
        CFfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', generictermreward, 'approxmethod', upperNoDisc}; % use this to not use KG* for terminal reward at time 'infinity'
        CFparamvec = { 't0', .1, 'tEND', THoriz, 'precfactor', 10, 'ceilfactor', 1.5, 'BaseFileName', [PDEnodiscbase PDEonbase] , 'matdir' , PDEmatfilebase, 'finiteT', finiteT  };
    else
        param.tEND = 1000 + param.t0;
        THoriz = param.tEND;
        warning('online learning requested for infinite horizon without discounting: not allowed -- try discounting or finite horizon: here, assuming finite horizon with Thoriz = %d',tHoriz);
        CFfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', generictermreward, 'approxmethod', upperNoDisc}; % use this to have KG* type rule at time 'infinity' for ca
%        CFfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', CFApproxValuefunc, 'approxmethod', upperNoDisc}; % use this to have KG* type rule at time 'infinity' for ca
        CFparamvec = { 't0', param.t0, 'tEND', THoriz, 'precfactor', 10, 'ceilfactor', 1.5, 'BaseFileName', [PDEnodiscbase PDEonbase] , 'matdir' , PDEmatfilebase, 'finiteT', finiteT  };
%        CFparamvec = { 't0', .1, 'tEND', 100000, 'precfactor', 10, 'ceilfactor', 1.5, 'BaseFileName', [PDEnodiscbase PDEonbase] , 'matdir' , PDEmatfilebase, 'finiteT', finiteT  };
        %figdir Figure\,  UnkVariance 0
    end
    scalevec = CFscalevec; 
    paramvec = [CFparamvec, CFfunctionset, baseparams];
    [scale, param] = PDEInputConstructor( scalevec, paramvec );
    [scale, param] = PDEInputValidator( scale, param )
    s0 = 1/(scale.gamma * param.t0)
    sEND = 1/(scale.gamma * param.tEND)
    tic
    [rval, MAXFiles] = PDESolnCompute(scale, param);
    % Load in the data structures form those computations
    toc

    % now turn to discounted version: online is ok for both finite and
    % infinite horizon
    % Set up generic functions for positive discounting
    CGApproxValuefunc=@(wvec,s,p1,p2) PDECGApproxValue(wvec,s,p1);   % this is valid terminal reward for discounted rewards, valued in time s currency
    upperDisc=@(s,p1,p2) CGApproxBoundW(s);
    CGscalevec = {'c', 0, 'sigma', 10e5, 'discrate', 0.0002, 'P', 1 };
    if finiteT
        CGfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', generictermreward, 'approxmethod', upperDisc};
        CGparamvec = { 't0', 0.002, 'tEND', THoriz, 'precfactor', 8, 'ceilfactor', 2.0, 'BaseFileName', [PDEdiscbase PDEonbase] , 'matdir' , PDEmatfilebase};
    else
        CGfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', CGApproxValuefunc, 'approxmethod', upperDisc};
        CGparamvec = { 't0', 0.002, 'tEND', 400000, 'precfactor', 8, 'ceilfactor', 2.0, 'BaseFileName', [PDEdiscbase PDEonbase] , 'matdir' , PDEmatfilebase};
        %figdir Figure\,  UnkVariance 0
    end

    scalevec = CGscalevec; 
    paramvec = [CGparamvec, CGfunctionset, baseparams];
    [scale, param] = PDEInputConstructor( scalevec, paramvec );
    [scale, param] = PDEInputValidator( scale, param )
    tic
    [rval, MAXFiles] = PDESolnCompute(scale, param);
    toc
        
end


end
