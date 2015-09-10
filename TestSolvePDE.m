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
generictermreward=@(wvec,s,p1,p2) PDEsimplereward(wvec);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
baseparams = { 'online', 0, 'retire', 0, 'DoPlot', 1 };

CFApproxValuefunc=@(wvec,s,p1,p2) PDECFApproxValue(wvec,s,p1);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
upperNoDisc=@(s,p1,p2) CFApproxBoundW(s);
CFfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', generictermreward, 'approxmethod', upperNoDisc}; % use this to not use KG* for terminal reward at time 'infinity'
CFfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', CFApproxValuefunc, 'approxmethod', upperNoDisc}; % use this to have KG* type rule at time 'infinity' for ca
CFscalevec = {'c', 1, 'sigma', 10e5, 'discrate', 0, 'P', 1};
CFparamvec = { 't0', .1, 'tEND', 100000, 'precfactor', 8, 'BaseFileName', 'CF' };
%figdir Figure\, matdir Matfiles\ UnkVariance 0

% Set up generic functions for positive discounting
CGApproxValuefunc=@(wvec,s,p1,p2) PDECGApproxValue(wvec,s,p1);   % this is valid terminal reward for discounted rewards, valued in time s currency
upperDisc=@(s,p1,p2) CGApproxBoundW(s);
CGfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', generictermreward, 'approxmethod', upperDisc};
CGfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', CGApproxValuefunc, 'approxmethod', upperDisc};
CGscalevec = {'c', 0, 'sigma', 10e5, 'discrate', 0.0002, 'P', 1 };
CGparamvec = { 't0', 0.002, 'tEND', 400000, 'precfactor', 6, 'BaseFileName', 'CG' };

% generic functions when the 'guesses' are still being made regarding the
% upper boundary's approximate value.
Guessdw=0.06; GuessNumW=1500; % these specific values are appropriate for case of P=1, c=1, discrate = 0.0002, for example
upperguessNoClue = [Guessdw GuessNumW]; % guesses for initial dw size, and for number of grid points above and below 0
Guessfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', generictermreward, 'approxmethod', upperguessNoClue};
Guessscalevec = {'c', 0, 'sigma', 10e5, 'discrate', 0.0001, 'P', 1 };
Guessparamvec = { 't0', 1, 'tEND', 20000, 'precfactor', 6, 'BaseFileName', 'Guess' };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                               %% 
%%       Create functions for use with zero discounting          %%   offline case
%%                                                               %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
BaseFileName = strcat(param.matdir,param.BaseFileName); % note, we wish to allow loading files by name without having the full solution or the full param: just the name and range of blocks to load
%[rval, cfSoln] = PDESolnLoad(BaseFileName,1,MAXFiles);
[rval, cfSoln] = PDESolnLoad(BaseFileName); % by default load all subgroups
if cfSoln.Header.PDEparam.DoPlot % do a bunch of diagnostics plots, save the eps files
    if ~exist('fignum','var'), fignum = 20; end;
    fignum = UtilPlotDiagnostics(fignum, cfSoln);
end

if ~exist('fignum','var'), fignum = 20; end;
[ rval, fignum, ~ ] = DoCFPlots( fignum, cfSoln );

%%% To TEST KGs stuff for case of no discounting
if ~exist('fignum','var'), fignum = 20; end;
s0 = 1/(scale.gamma*param.tEND);
dw = .0005;
bigw = 40;
wvec = (-bigw:bigw)*dw;
[voikg, dsvec]=PDECFApproxValue(wvec,s0,scale);
voikg = voikg - max(0,wvec);
fignum=fignum+1;figure(fignum);
plot(wvec,voikg);
fignum=fignum+1;figure(fignum);
plot(wvec/scale.beta,1./(scale.beta*dsvec));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                               %% 
%%     Create functions for use with positive discounting        %%  offline case
%%                                                               %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scalevec = CGscalevec; 
paramvec = [CGparamvec, CGfunctionset, baseparams];
[scale, param] = PDEInputConstructor( scalevec, paramvec );
[scale, param] = PDEInputValidator( scale, param )
tic
[rval, MAXFiles] = PDESolnCompute(scale, param);
toc
% Load in the data structures form those computations
BaseFileName = strcat(param.matdir,param.BaseFileName); % note, we wish to allow loading files by name without having the full solution or the full param: just the name and range of blocks to load
%[rval, cgSoln] = PDESolnLoad(BaseFileName,1,MAXFiles);
[rval, cgSoln] = PDESolnLoad(BaseFileName);
if cgSoln.Header.PDEparam.DoPlot % do a bunch of diagnostics plots, save the eps files
    if ~exist('fignum','var'), fignum = 20; end;
    fignum = UtilPlotDiagnostics(fignum, cgSoln);
end

if ~exist('fignum','var'), fignum = 20; end;
[ rval, fignum, ~ ] = DoCGPlots( fignum, cgSoln ); % can pass with only one argument, in which case Matfiles\CG0.mat is checked for loading in pde solution
% alternatively, DoCGPlots can take a string, or can be left blank to load
% in files in the default location
%[ rval, figout, pdeSolnStruct ] = DoCGPlots( fignum, 'Matfiles\CG' );
%[ rval, figout, pdeSolnStruct ] = DoCGPlots( fignum );

%%% To TEST KGs stuff for case of discounting
if ~exist('fignum','var'), fignum = 20; end;
s0 = 1/scale.gamma/param.tEND;
dw = .01;
bigw = 2000;
wvec = (-2*bigw:bigw)*dw;
[voikg, dsvec]=PDECGApproxValue(wvec,s0,scale);     % with 'scale' passed, a number of powers of 2 of budgets are checked: might be better to do without 'scale' passed
voikg = voikg - max(0,wvec);
[voikg2, dsvec2]=PDECGApproxValue(wvec,s0);           % without 'scale' passed, a large number of values of 's' are checked.
voikg2 = voikg2 - max(0,wvec);
fignum=fignum+1;figure(fignum);
plot(wvec,voikg,'-',wvec,voikg2,'-.');
legend('with scale', 'without');
fignum=fignum+1;figure(fignum);
%plot(wvec/scale.beta,1./(scale.gamma*dsvec));
%plot(wvec/scale.beta,1./(scale.gamma*(s0-dsvec)) - 1./(scale.gamma*s0));
plot(wvec,dsvec,'-',wvec,dsvec2,'-.');
legend('with scale', 'without');

[kg, dsvec]=PDECGKGs(wvec,s0,scale);     % with 'scale' passed, a number of powers of 2 of budgets are checked: might be better to do without 'scale' passed
[kg2, dsvec2]=PDECGKGs(wvec,s0);           % without 'scale' passed, a large number of values of 's' are checked.
fignum=fignum+1;figure(fignum);
plot(wvec,kg,'-',wvec,kg2,'-.');
legend('with scale', 'without');
fignum=fignum+1;figure(fignum);
%plot(wvec/scale.beta,1./(scale.gamma*dsvec));
plot(wvec,dsvec,'-',wvec,dsvec2,'-.');
legend('with scale', 'without');



t0 = param.t0
s0 = 1/scale.gamma/t0
dw = .1;
bigw = 10;
wvec = (-2*bigw:bigw)*dw/scale.beta;
lookahead = 1;
kgfrac = wvec.^2 * t0 / scale.sigma^2;
kgs = max( 0, (t0/4) * (kgfrac - 1 + sqrt(kgfrac.^2 + 6*kgfrac + 1)));
lookahead = kgs;
zvar = scale.sigma^2 * lookahead / (t0 * (t0+lookahead));

maxv=-scale.c*lookahead + sqrt(zvar) .* PsiNorm(-wvec ./ sqrt(zvar));


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
    fignum = UtilPlotDiagnostics(fignum,GuessSoln);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD STUFF %%%%%%%%%%%%%%%%%%%%%%
m=0;         % value of retirement option in (y,t) space
mu0=0;

wvec

APrioriRegret=(scale.sigma/sqrt(param.t0)) * PsiNorm( abs(mu0-m) / (scale.sigma/sqrt(param.t0)))
tmppp=PDEGetVals(cfSoln,(mu0-m)*scale.beta,1/(param.t0*scale.gamma))/scale.beta
APrioriRegret-tmppp


