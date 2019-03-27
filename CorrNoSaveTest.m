% SolvePlotCFCG.m:
% A macro with the code required to:
%   1. Set up the path structure in matlab which allows the general PDE
%   solution code to be called.
%   2. Create files with the solutions to the free boundary PDE problems in
%   Chick & Frazier and in Chick & Gans, for offline learning, and an
%   extension which accounts for online learning (requires c=1 for C&F and
%   requires additional conditions for C&G, to be able to use the online
%   learning solution).
%   3. Loads the files with the solutions into memory, and generates a
%   number of diagnostic plots for those solutions.
%   4. Generates the plots for the C&F paper and for the C&G paper.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                               %% 
%%     SETUP: For files and directory, to be called every time.
%%                                                               %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PDELocalSetPaths(); % set up paths to PDE package
PDELocalInit;       % set up some variables (all starting with PDE - avoid variables with such names
[cgSoln, cfSoln, cgOn, cfOn] = PDELoadSolnFiles(PDEmatfilebase, onflag);

fignum = 20;

onflag = true;
onflag = false;
THoriz = 1000;    % set to below 0 for infinite horizon. Look ahead this many steps for adaptive stopping time.

if THoriz <= 0;
    finiteT = false; % by default, compute infinite horizon solution
else
    finiteT = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                               %% 
%%     SETUP Initial setup stuff with zero discount, independent rewards
%%                                                               %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up generic functions for zero discounting,
% these functions should return two vectors of the same size as wvec: the
% first should be the expected reward at time s given wvec, the second
% should be the expected number of samples to take additionally to achieve
% it.
% First set things up for offline learning, then do online learning files
% if needed.

generictermreward=@(wvec,s,p1,p2) PDEsimplereward(wvec);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
CFApproxValuefunc=@(wvec,s,p1,p2) PDECFApproxValue(wvec,s,p1);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
upperNoDisc=@(s,p1,p2) CFApproxBoundW(s);

baseparams = { 'online', 0, 'DoPlot', 1, 'DoFileSave', 0 }; % NOTE: DoSaveFile is default by true, 

myk = 4;  % number of alternatives
myt0 = 5; %


c = 1;          % cost per sample
mysig = 10e5;   % sample variance
P = 1;      % population size upon finishing
discrate = 0.0;

myi = 2;        % set this to be the alternative to be tested
mytt = myt0;    % set this to be the effectivenumber of samples of alternative myi
myMutVec = 1:myk; % this should be the vector of means given samples so far
mySigtMat = mysig^2 * (eye(myk) + ones(myk)/10);
mySigtMat = mysig^2 * cov(randn(50*myk,myk));
eveci = 0*myMutVec; eveci(myi)=1;

CFscalevec = {'c', c, 'sigma', mysig, 'discrate', discrate, 'P', P};
distribparams = { 'myi', myi, 'mytt', mytt, 'myMutVec', myMutVec, 'mySigtMat', mySigtMat }; 

%if finiteT
%    CFfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', generictermreward, 'approxmethod', upperNoDisc}; % use this to not use KG* for terminal reward at time 'infinity'
%    CFparamvec = { 't0', myt0, 'tEND', myt0+THoriz, 'precfactor', 10, 'ceilfactor', 1.1, 'BaseFileName', [PDEnodiscbase PDEoffbase] , 'matdir' , PDEmatfilebase, 'finiteT', finiteT  };
%else
    CFfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', CFApproxValuefunc, 'approxmethod', upperNoDisc}; % use this to have KG* type rule at time 'infinity' for ca
%    CFparamvec = { 't0', .1, 'tEND', 100000, 'precfactor', 10, 'ceilfactor', 1.1, 'BaseFileName', [PDEnodiscbase PDEoffbase] , 'matdir' , PDEmatfilebase, 'finiteT', finiteT  };
    CFparamvec = { 't0', myt0, 'tEND', myt0+THoriz, 'precfactor', 12, 'ceilfactor', 1.15, 'BaseFileName', [PDEnodiscbase PDEoffbase] , 'matdir' , PDEmatfilebase, 'finiteT', finiteT  };
    %figdir Figure\,  UnkVariance 0
%end

scalevec = CFscalevec;
paramvec = [CFparamvec, CFfunctionset, baseparams, distribparams];
[scale, param] = PDEInputConstructor( scalevec, paramvec );
[scale, param] = PDEInputValidator( scale, param )
s0 = 1/(scale.gamma * param.t0)
sEND = 1/(scale.gamma * param.tEND)
tic
%[rval, MAXFiles] = PDESolnCompute(scale, param);
[rval, output] = PDESolnCompute(scale, param);
% Load in the data structures form those computations
toc

mu0stotest = (-3:0.1:3)*mysig/sqrt(param.t0);
Bwsval = (interp2(output.svec,output.wvec,output.Bwsmatrix,s0,scale.beta*mu0stotest)/scale.beta+max(mu0stotest',0))

if param.online 
    valfromPDE = PDEGetVals(cfOn,scale.beta*mu0stotest,s0)/scale.beta;
else
    valfromPDE = PDEGetVals(cfSoln,scale.beta*mu0stotest,s0)/scale.beta;
end

diffInfHoriztoFiniteT = valfromPDE - Bwsval;
rationInfHoriztoFiniteT = (valfromPDE - Bwsval)./valfromPDE;

VvecKG = param.approxvaluefunc(scale.beta*mu0stotest,s0,cfSoln.Header.PDEscale,cfSoln.Header.PDEparam)/scale.beta;

figure(fignum); fignum=fignum+1;
plot(mu0stotest,valfromPDE,'--',mu0stotest,Bwsval,'-.',mu0stotest,max(0,mu0stotest),'-',mu0stotest,VvecKG,'.')
legend('PDE from file','quick comp pde','hockey stick','KGval')


