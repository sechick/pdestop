function [CGOff, CFOff, CGOn, CFOn] = PDECreateSolnFiles(fName, onflag, THoriz)
% PDECreateSolnFiles: Solves free boundary problem in scaled coordinates for the
% sequential sampling for simulation optimization. 
%
% For sample calling conventions, including how to set up the fields of the
% structures PDEscale and PDEparam, please see TestSolvePDE.m.
%
% INPUTS: 
%   fName: base name for files to hold the data structures with the
%       solutions of the free boundary PDE.
%   onflag: if true, then online learning solutions are also created, if
%       false, then only offline learning solutions are created.
%   THoriz: set to positive, finite value if finite time horizon is desired, 
%       with given number of replications. If 0 or negative, then infinite 
%       time horizon is assumed. THoriz set assuming params scaled so that
%       sigma = 1, offline learning, and (c=0, disc=10^-3) for discounted
%       and (c=1,disc=0) for undiscounted.
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
%   Side effect: files of name fName[On/Off]<n>.mat for n=1,2,...,numfiles
%       which contain solutions for PDE in standardized coordinates 
%       Contents include: upper and lower boundaries, standardized: reward
%       function, PCS, expected number of samples to stop time, ...
%       also creates fName[On/Off]0.mat, with summary information for files.
%           base name (with On/Off) can be found in <structurename>.Header.fName
%   Side effect: a bunch of plots as calculations are run - if they don't
%       'look right' then there is probably a numerical stability issue to
%       address.
%
% Online learning: not yet supported.
% CFP: not yet supported.
%
% Code provided 'as is' with no warrantees of correctness.
%
% Author: S Chick, 2015, all rights reserved.
% Created: 2015 08 17.

if nargin < 1
    fName = 'PDE-';
end
if nargin < 2
    onflag = false;  % default is to only generate files for offline learning
end
if nargin < 3
    THoriz = -1;    % default is to run infinite horizon analysis
end

sigma=1000;     % variance per single sample in (y,t) space
m=0.0;         % value of retirement option in (y,t) space
mu0=0.0;
% ** Terminal reward functions: What is expected value of stopping at time
% s with means in wvec
basicstopfunc=@(wvec,s,p1,p2) max(wvec,0);   % this is valid terminal reward for discounted rewards, valued in time s currency
% ** Approximate value to go function: 
CFApproxValuefunc=@(wvec,s,p1,p2) PDECFKGs(wvec,s);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
CGApproxValuefunc=@(wvec,s,p1,p2) PDECGKGs(wvec,s);   % this is valid terminal reward for discounted rewards, valued in time s currency
upperNoDisc=@(s,p1,p2) CFApproxBoundW(s);
upperDisc=@(s,p1,p2) CGApproxBoundW(s);

PDEparam.online = false;% true if online learning is active, false otherwise
PDEparam.retire = m;    % value of retirement option: ignored for discrate=0, taken to be alternative for discrate>0

if THoriz <= 0;
    PDEparam.finiteT = false; % by default, compute infinite horizon solution
else
    PDEparam.finiteT = true;
end

fileNameModifiers = {'Off', 'On'};
if onflag
    loopmax = 2;
else % if online offline learning needed, set structure for online learning to empty structure
    loopmax = 1;
    CGOn=[];
    CFOn=[];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:loopmax
    currentSuffix = fileNameModifiers(i);
    baseName = sprintf('%s%s',fName,currentSuffix{:});
    if PDEparam.finiteT 
        baseName = [baseName 'Fin'];
    end
    PDEparam.online = (i > 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % first do the case of zero discount rate
    c=1;         % cost per single sample in (y,t) space
    discrate = 0.00;% zero discount rate per sample for C&F paper (use >0 for C&G paper)
    t0=0.01;%0.25; %.25
    if THoriz > 0
        tEND = max(2*t0, THoriz);   % make sure there is a nonempty range
    else
        tEND = 1000;  %5000
    end
        
    PDEscale.c=c;           % cost per single sample in (y,t) space
    PDEscale.sigma=sigma;	% standard deviation of a single sample in (y,t) space
    PDEscale.P=1;           % multiplier times the mean reward for the adoption decision
    PDEscale.discrate = discrate;% zero discount rate per sample for C&F paper (use >0 for C&G paper)
    PDEscale = PDEScaleStandardize(PDEscale);    % compute alpha, beta, gamma, kappainv for this problem

   
    PDEparam.BaseFileName = ['CF' baseName];      % lower range of values of t0 to compute, where unknown mean W ~ Normal(mu0, sigma^2/t0)
    PDEparam.t0 = t0;      % lower range of values of t0 to compute, where unknown mean W ~ Normal(mu0, sigma^2/t0)
    PDEparam.tEND = tEND;  % upper range of values of t0 to compute, where unknown mean W ~ Normal(mu0, sigma^2/t0)
    PDEparam.termrewardfunc = basicstopfunc;
    PDEparam.approxvaluefunc = CFApproxValuefunc ;
    PDEparam.approxmethod = upperNoDisc;
    PDEparam.precfactor = 10.0;	% increase to increase the precision of the 'dw' grid in normed coordinates. must be at least 1.0.
    PDEparam.DoPlot = true;       % true to turn on diagnostic plots

    tic
    [rval, MAXFiles] = PDESolnCompute(PDEscale, PDEparam);
    if i==1
        [rval, CFOff] = PDESolnLoad(BaseNameFile,1,MAXFiles);
    else
        [rval, CFOn] = PDESolnLoad(BaseNameFile,1,MAXFiles);
    end
    toc
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now do the [positive discount rate
    c=0;
    discrate = 0.001;% zero discount rate per sample for C&F paper (use >0 for C&G paper)
    t0=1;%0.25; %.25
    if THoriz > 0
        tEND = max(2*t0, THoriz);   % make sure there is a nonempty range
    else
        tEND = 1000;  %5000
    end
    
    PDEscale.c=c;           % cost per single sample in (y,t) space
    PDEscale.sigma=sigma;	% standard deviation of a single sample in (y,t) space
    PDEscale.P=1;           % multiplier times the mean reward for the adoption decision
    PDEscale.discrate = discrate;% zero discount rate per sample for C&F paper (use >0 for C&G paper)
    PDEscale = PDEScaleStandardize(PDEscale);    % compute alpha, beta, gamma, kappainv for this problem

    BaseNameFile=['CG' baseName];
    PDEparam.termrewardfunc = basicstopfunc;
    PDEparam.approxvaluefunc = CGApproxValuefunc;
    PDEparam.approxmethod = upperDisc;
    PDEparam.precfactor = 6.0;	% increase to increase the precision of the 'dw' grid in normed coordinates. must be at least 1.0.
    PDEparam.DoPlot = true;       % true to turn on diagnostic plots

    tic
    [rval, MAXFiles] = PDESolnCompute(BaseNameFile, PDEscale, PDEparam);
    if i==1
        [rval, CGOff] = PDESolnLoad(BaseNameFile,1,MAXFiles);
    else
        [rval, CGOn] = PDESolnLoad(BaseNameFile,1,MAXFiles);
    end
    toc

end
