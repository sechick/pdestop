%Test recursion in (w,s) space for sup_S E[ max(W_S/S,0) - 1/S | w_0, s_0]
%problem (the
%normalized problem of maximizing value of decision (pick best of 0 or
%posterior reward) plus sampling cost.  In (y,t) space, we have:
%   sup_T[ max(Y_T/T,0) - c T | y_0/t_0 = \mu_0, t_0] + c t_0,
%for nonanticipating stopping rules T>t_0, and where t = t_0 + n, n is
%number of samples, and y_t = y_0 + \sum_{i=1}^n x_i, X_i ~
%Normal(A,\sigma^2), and prior on A ~ Normal(\mu_0, \sigma_0^2), and
%further t_0 = \sigma^2 / \sigma_0^2.

% Set up generic functions for zero discounting,
% these functions should return two vectors of the same size as wvec: the
% first should be the expected reward at time s given wvec, the second
% should be the expected number of samples to take additionally to achieve
% it.

PDELocalInit;

DelayOfflineReward=@(wvec,s,p1,p2) DelayOfflineSimpleReward(wvec,s,p1,p2);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
DelayOffApproxV=@(wvec,s,p1,p2) DelayNodiscOffApproxV(wvec,s,p1,p2);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
upperNoDisc=@(s,p1,p2) CFApproxBoundW(s);
DelayNodiscOffFuncset = {'termrewardfunc', DelayOfflineReward, 'approxvaluefunc', DelayOffApproxV, 'approxmethod', upperNoDisc}; % use this to have KG* type rule at time 'infinity' for ca
DelayNodiscOffFuncset = {'termrewardfunc', DelayOfflineReward, 'approxvaluefunc', DelayOfflineReward, 'approxmethod', upperNoDisc}; % use this to not use KG* for terminal reward at time 'infinity'

% Put in data for solution for Stent example
annualratepatients = 907;
tauval = 907;
Stentt0 = 20;
StentMaxSamps = 2000;
annualdiscrate = .10;
perpatientdiscrate = annualdiscrate / annualratepatients;
adoptionsize = 2e6;
samplecost = 200;
samplestdev = 17538;

% Try first when there is zero discount rate
StentScale = {'c', samplecost, 'sigma', samplestdev, 'discrate', 0, 'P', adoptionsize, 'tau', tauval };
StentParam = { 't0', Stentt0, 'tEND', Stentt0+StentMaxSamps, 'precfactor', 6, 'BaseFileName', sprintf('StentOffTau%d',tauval) };
baseparams = { 'online', 0, 'retire', 0, 'DoPlot', 1, 'finiteT', true };
%Force to be a finite time process with given time horizon, 
scalevec = StentScale; 
paramvec = [StentParam, DelayNodiscOffFuncset, baseparams];
[scale, param] = PDEInputConstructor( scalevec, paramvec );
[scale, param] = PDEInputValidator( scale, param )
tic
[rval, MAXFiles] = PDESolnCompute(scale, param);
% Load in the data structures form those computations
toc
BaseFileName = strcat(param.matdir,param.BaseFileName); % note, we wish to allow loading files by name without having the full solution or the full param: just the name and range of blocks to load
[rval, cfSoln] = PDESolnLoad(BaseFileName,1,MAXFiles);
if cfSoln.Header.PDEparam.DoPlot % do a bunch of diagnostics plots, save the eps files
	if ~exist('fignum','var'), fignum = 20; end;
    fignum = UtilPlotDiagnostics(fignum,cfSoln);
end



dw = .005;
bigw = 15;
wvec = (-bigw:bigw)*dw;
[voikg, dsvec]= cfSoln.Header.PDEparam.approxvaluefunc(wvec,s0,scale,param);
voikg = scale.P * (voikg - max(0,wvec));
fignum=fignum+1;figure(fignum);
plot(wvec,voikg);
fignum=fignum+1;figure(fignum);
plot(wvec/scale.beta,scale.P./(scale.gamma*dsvec));

m=0;         % value of retirement option in (y,t) space
mu0=0;

APrioriRegret=(scale.sigma/sqrt(param.t0)) * PsiNorm( abs(mu0-m) / (scale.sigma/sqrt(param.t0)))
tmppp=PDEGetVals(cfSoln,(mu0-m)*scale.beta,1/(param.t0*scale.gamma))/scale.beta
APrioriRegret-tmppp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                               %% 
%%     Create functions for use with positive discounting        %%  offline case
%%                                                               %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DelayOfflineReward=@(wvec,s,p1,p2) DelayOfflineSimpleReward(wvec,s,p1,p2);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
DelayDiscOffApproxValuefunc=@(wvec,s,p1,p2) DelayDiscOffApproxV(wvec,s,p1);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
upperDisc=@(s,p1,p2) CFApproxBoundW(s);
DelayDiscOffFuncset = {'termrewardfunc', DelayOfflineReward, 'approxvaluefunc', DelayDiscOffApproxValuefunc, 'approxmethod', upperDisc}; % use this to have KG* type rule at time 'infinity' for ca
DelayDiscOffFuncset = {'termrewardfunc', DelayOfflineReward, 'approxvaluefunc', DelayOfflineReward, 'approxmethod', upperDisc}; % use this to not use KG* for terminal reward at time 'infinity'
tauval = 907;
DStentScale = {'c', 200, 'sigma', 17538, 'discrate', 0.01, 'P', 2e6, 'tau', tauval };
DStentParam = { 't0', 20, 'tEND', 2000, 'precfactor', 6, 'BaseFileName', sprintf('DelDisOffTau%d',tauval) };
baseparams = { 'online', 0, 'retire', 0, 'DoPlot', 1, 'finiteT', true };
%Force to be a finite time process with given time horizon, 

scalevec = DStentScale; 
paramvec = [DStentParam, DelayDiscOffFuncset, baseparams];
[scale, param] = PDEInputConstructor( scalevec, paramvec );
[scale, param] = PDEInputValidator( scale, param )
tic
[rval, MAXFiles] = PDESolnCompute(scale, param);
toc
% Load in the data structures form those computations
BaseFileName = strcat(param.matdir,param.BaseFileName); % note, we wish to allow loading files by name without having the full solution or the full param: just the name and range of blocks to load
[rval, cgdSoln] = PDESolnLoad(BaseFileName,1,MAXFiles);
if cgdSoln.Header.PDEparam.DoPlot % do a bunch of diagnostics plots, save the eps files
    if ~exist('fignum','var'), fignum = 20; end;
    fignum = UtilPlotDiagnostics(fignum, cgdSoln);
end


%%% To TEST KGs stuff for case of no discounting
if ~exist('fignum','var'), fignum = 20; end;
s0 = 1/(scale.gamma*param.t0)
sEND = 1/(scale.gamma*param.tEND)
Icost = 0;
[rval, cgSoln] = PDESolnLoad([PDEmatfilebase PDEdiscbase],1,6);
[ Vvec, Bwsval, ENvec, PCSvec ] = PDEGetVals(cgSoln,wvec - beta*Icost / scale.P + beta*scale.c / (scale.P * scale.discrate),s0) 
[ Vvec, Bwsval, ENvec, PCSvec ] = PDEGetVals(cgdSoln,wvec - beta*Icost / scale.P + beta*scale.c / (scale.P * scale.discrate),s0) 
Vvec = scale.P * Vvec / beta - scale.c/scale.discrate;
muvec = wvec / beta;



%%% To TEST KGs stuff for case of discounting
if ~exist('fignum','var'), fignum = 20; end;
s0 = 1/scale.gamma/param.tEND;
dw = .01;
bigw = 20;
wvec = (-2*bigw:bigw)*dw;
[voikg, dsvec]=PDECGApproxValue(wvec,s0,scale);
voikg = voikg - -max(0,wvec);
plot(wvec,voikg);
fignum=fignum+1;figure(fignum);
plot(wvec,voikg);
fignum=fignum+1;figure(fignum);
plot(wvec/scale.beta,1./(scale.gamma*dsvec));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                               %% 
%%       Create functions for use with ad hoc grid size and discounting          %%   offline case
%%                                                               %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scalevec = Guessscalevec; 
paramvec = [Guessparamvec, Guessfunctionset, baseparams];
[scale, param] = PDEInputConstructor( scalevec, paramvec );
[scale, param] = PDEInputValidator( scale, param )
tic
[rval, MAXFiles] = PDESolnCompute(scale, param);
toc
% Load in the data structures form those computations
BaseFileName = strcat(param.matdir,param.BaseFileName); % note, we wish to allow loading files by name without having the full solution or the full param: just the name and range of blocks to load
[rval, GuessSoln] = PDESolnLoad(BaseFileName,1,MAXFiles);
if GuessSoln.Header.PDEparam.DoPlot % do a bunch of diagnostics plots, save the eps files
     if ~exist('fignum','var'), fignum = 20; end;
     fignum = UtilPlotDiagnostics(fignum, GuessSoln);
end




