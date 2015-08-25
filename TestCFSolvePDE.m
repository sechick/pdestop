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
CFApproxValuefunc=@(wvec,s,p1,p2) PDECFKGs(wvec,s,p1);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
upperNoDisc=@(s,p1,p2) CFApproxBoundW(s);
CFfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', CFApproxValuefunc, 'approxmethod', upperNoDisc}; % use this to have KG* type rule at time 'infinity' for ca
CFfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', generictermreward, 'approxmethod', upperNoDisc}; % use this to not use KG* for terminal reward at time 'infinity'
CFscalevec = {'c', 1, 'sigma', 10e5, 'discrate', 0, 'P', 1};
CFparamvec = { 't0', .3, 'tEND', 20000, 'precfactor', 6, 'BaseFileName', 'CF' };
baseparams = { 'online', 0, 'retire', 0, 'DoPlot', 1 };
%figdir Figure\, matdir Matfiles\ UnkVariance 0

% Set up generic functions for positive discounting
CGApproxValuefunc=@(wvec,s,p1,p2) PDECGKGs(wvec,s,p1);   % this is valid terminal reward for discounted rewards, valued in time s currency
upperDisc=@(s,p1,p2) CGApproxBoundW(s);
CGfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', generictermreward, 'approxmethod', upperDisc};
CGfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', CGApproxValuefunc, 'approxmethod', upperDisc};
CGscalevec = {'c', 0, 'sigma', 10e5, 'discrate', 0.0002, 'P', 1 };
CGparamvec = { 't0', 1, 'tEND', 20000, 'precfactor', 6, 'BaseFileName', 'CG' };

% generic functions when the 'guesses' are still being made regarding the
% upper boundary's approximate value.
Guessdw=0.06; GuessNumW=1500; % these specific values are appropriate for case of P=1, c=1, discrate = 0.0002, for example
upperguessNoClue = [Guessdw GuessNumW]; % guesses for initial dw size, and for number of grid points above and below 0
Guessfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', generictermreward, 'approxmethod', upperguessNoClue};
Guessscalevec = {'c', 0, 'sigma', 10e5, 'discrate', 0.0002, 'P', 1 };
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
tic
[rval, MAXFiles] = PDESolnCompute(scale, param);
% Load in the data structures form those computations
toc
BaseFileName = strcat(param.matdir,param.BaseFileName); % note, we wish to allow loading files by name without having the full solution or the full param: just the name and range of blocks to load
[rval, cfSoln] = PDESolnLoad(BaseFileName,1,MAXFiles);
if cfSoln.Header.PDEparam.DoPlot % do a bunch of diagnostics plots, save the eps files
    UtilPlotDiagnostics(cfSoln);
end

%%% To TEST KGs stuff for case of no discounting
if ~exist('fignum','var'), fignum = 20; end;
s0 = 1/(scale.gamma*param.tEND);
dw = .005;
bigw = 15;
wvec = (-bigw:bigw)*dw;
[voikg, dsvec]=PDECFKGs(wvec,s0,scale);
voikg = voikg - -max(0,wvec);
fignum=fignum+1;figure(fignum);
plot(wvec,voikg);
fignum=fignum+1;figure(fignum);
plot(wvec/scale.beta,1./(scale.gamma*dsvec));

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
[rval, cgSoln] = PDESolnLoad(BaseFileName,1,MAXFiles);
if cfSoln.Header.PDEparam.DoPlot % do a bunch of diagnostics plots, save the eps files
    UtilPlotDiagnostics(cgSoln);
end

%%% To TEST KGs stuff for case of discounting
if ~exist('fignum','var'), fignum = 20; end;
s0 = 1/scale.gamma/param.tEND;
dw = .01;
bigw = 20;
wvec = (-2*bigw:bigw)*dw;
[voikg, dsvec]=PDECGKGs(wvec,s0,scale);
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
    UtilPlotDiagnostics(GuessSoln);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD STUFF %%%%%%%%%%%%%%%%%%%%%%
m=0;         % value of retirement option in (y,t) space
mu0=0;

APrioriRegret=(scale.sigma/sqrt(param.t0)) * PsiNorm( abs(mu0-m) / (scale.sigma/sqrt(param.t0)))
tmppp=PDEGetVals(cfSoln,(mu0-m)*scale.beta,1/(param.t0*scale.gamma))/scale.beta
APrioriRegret-tmppp

% start to reconstruct the values for plotting....
sbvec = cfSoln.Computed.accumsvec;
up = cfSoln.Computed.accumupper;
lo = cfSoln.Computed.accumlower;
findex = cfSoln.Computed.fileindx;

if param.DoPlot
    if ~exist('fignum','var'), fignum = 20; end;
    if isa(cfSoln.Header.PDEparam.approxmethod,  'function_handle')
    %    dt = 1/scale.gamma/s0 - 1/scale.gamma/(s0+ds)
    %    dmu=dw/beta
        mysmallfontsize = 14;
        myfontsize = 16;
        fignum=fignum+1;figure(fignum);
        tstcurv=cfSoln.Header.PDEparam.approxmethod(sbvec);
        plot(sbvec,up,'--',sbvec,tstcurv,'-.');set(gca,'FontSize',mysmallfontsize);
        legend('From PDE','Quick Approx.','Location','SouthEast')
        xlabel('Reverse time scale,  s=1/(\gamma t)','FontSize',myfontsize,'FontName','Times'); ylabel('Rescaled mean, w','FontSize',myfontsize,'FontName','Times')
        %title('Upper stopping boundary','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat('BoundaryWS','.eps');
        print('-deps',mytitle);	

        betastepvec=[1 5 10 20 50 100 200 500 1000];
        numbetas=length(betastepvec);
        betamatrix=zeros(numbetas,length(up));

        fignum=fignum+1;figure(fignum);
        loglog(1/scale.gamma./sbvec,up/scale.beta,'--',1/scale.gamma./sbvec,tstcurv/scale.beta,'-.');set(gca,'FontSize',mysmallfontsize);
        legend('From PDE','Quick Approx.','Location','NorthEast')
        xlabel('Effective number of replications, n_t','FontSize',myfontsize,'FontName','Times'); ylabel('Posterior mean, y_t/n_t','FontSize',myfontsize,'FontName','Times')
        %title('Upper stopping boundary','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat('BoundaryYT','.eps');
        print('-deps',mytitle);	
    end
end


    mwig = param.retire*scale.beta;
    
    % GET INDEX INTO RIGHT STUCTURE BASED ON TIME
    [a, b] = min(s0 > sbvec);
    ind=findex(b)
    wvec = cfSoln.Data(ind).wvec;
    svec = cfSoln.Data(ind).svec;
    Bwsmatrix = cfSoln.Data(ind).Bwsmatrix;
    Bmu0t0=mwig/scale.beta+interp2(svec,wvec,Bwsmatrix,s0,w0)/scale.beta
    Bw0s0=mwig + interp2(svec,wvec,Bwsmatrix,s0,w0)
    % This is the VOI info - need to add in the stopping to get the full
    % value function
end

PDEGetVals(cfSoln,wvec,sval) % try to get solution for values of w in wvec at time sval, given pde solution cfSoln



%%%%%%%%% convert to standardized version in (w,s) with w brownian motion
%%%%%%%%% in -s time scale
mwig = scale.beta * param.retire
w0= scale.beta*(mu0-param.retire)
s0 = 1 / (scale.gamma * param.t0)
sEND = 1 / (scale.gamma * param.tEND)

% set up info for plotting stuff
myfontsize=16
mysmallfontsize=14
points = 144*3 %spacing between labels on contours - made so that only one label appears per line


%[ta tb tc]=fminbnd('Bztauapproxc',0,1/scale.gamma/n0base,optimset('TolFun',scale.beta/scale.sigma,'MaxIter',10^4),0,1/scale.gamma/n0base);       % FIXED VERSION
[ta tb tc]=fminbnd('Bztauapproxc',0,s0,optimset('TolFun',scale.beta/scale.sigma,'MaxIter',10^4),0,s0);       % FIXED VERSION
-tb/scale.beta - scale.c/scale.gamma
sig0 = scale.sigma / sqrt(param.t0)
sig0*PsiNorm(-mu0/sig0) - (-tb/scale.beta - scale.c/scale.gamma)

%ds
%dt = 1/scale.gamma/s0 - 1/scale.gamma/(s0+ds)
%dw
%dmu=dw/scale.beta

%Bmu0t0=m+interp2(svec,wvec,Bwsmatrix,s0,w0)/scale.beta
%Bw0s0=scale.beta*m + interp2(svec,wvec,Bwsmatrix,s0,w0)

% CHUNK5: Generate some plots
% Convert from (w,s) coordinates to (y/t, t) coordinates (NOT (y,t)
% coordinates!!!!)

betastepvec=[1 5 10 20 50 100 200 500 1000];
numbetas=length(betastepvec);
betamatrix=zeros(numbetas,length(up));

if ~exist('fignum','var'), fignum = 20; end;
fignum=fignum+1;figure(fignum);
tstcurv=cfSoln.Header.PDEparam.approxmethod(sbvec);
loglog(1/scale.gamma./sbvec,up/scale.beta,'--',1/scale.gamma./sbvec,tstcurv/scale.beta,'-.');set(gca,'FontSize',mysmallfontsize);
legend('From PDE','Quick Approx.','Location','NorthEast')
xlabel('Effective number of replications, n_t','FontSize',myfontsize,'FontName','Times'); ylabel('Posterior mean, y_t/n_t','FontSize',myfontsize,'FontName','Times')
%title('Upper stopping boundary','FontSize',myfontsize,'FontName','Times')
mytitle = strcat('BoundaryYT','.eps');
print('-deps',mytitle);	

if ~exist('fignum','var'), fignum = 20; end;
fignum=fignum+1;figure(fignum);
tstcurv= cfSoln.Header.PDEparam.approxmethod(sbvec);
xlabel('s coord');
title('ratio of upper bound from grid divided by approx to bound');
plot(sbvec,up./tstcurv)

% Try to get one-step estimate of boundary
for j=1:numbetas;
    ybndup1step = up/scale.beta; %initialize vector
    tvec = 1/scale.gamma./sbvec;
    ytinit = ybndup1step(length(ybndup1step))*tvec(length(ybndup1step));
    betamatrix(j,:) = ytinit;
    nrepslookahead=betastepvec(j)
    for i=length(ybndup1step):-1:1
        [ytinit, fval, exitflag]=fzero(@PsiNormRepsRoot,ytinit,optimset('TolX',1e-8),scale.sigma,tvec(i),nrepslookahead,nrepslookahead*scale.c);
        ybndup1step(i)=ytinit;
        ytinit=ytinit*0.99;
    end
    betamatrix(j,:)=ybndup1step;
end
kgstar=max(betamatrix);

if ~exist('fignum','var'), fignum = 20; end;
fignum=fignum+1;figure(fignum);
%loglog(1/scale.gamma./sbvec,up/scale.beta,'--',1/scale.gamma./sbvec,tstcurv/scale.beta,'-.',tvec,ybndup1step./tvec,':o',tvec,ybndupNstep./tvec,'-x')
%legend('From PDE','Quick Approx.','One-step lookahead','N-step lookahead','Location','NorthEast')
tstcurv=cfSoln.Header.PDEparam.approxmethod(sbvec);
loglog(1/scale.gamma./sbvec,tstcurv/scale.beta,'-',1/scale.gamma./sbvec,up/scale.beta,'-',1/scale.gamma./sbvec,kgstar./tvec,'--',tvec,betamatrix(1,:)./tvec,':',tvec,betamatrix(3,:)./tvec,'.',tvec,betamatrix(6,:)./tvec,'-.')
legend('Quick approx','PDE','KG_*','KG_1 = 1-step lookahead','KG_{10} = 10-step lookahead','KG_{100} =100-step lookahead','Location','SouthWest')
%loglog(1/scale.gamma./sbvec,tstcurv/scale.beta,'-.',tvec,ybndup1step./tvec,':')
%legend('PDE or Quick Approx.','One-step lookahead','Location','NorthEast')
xlabel('Effective number of replications, n_t','FontSize',myfontsize,'FontName','Times'); ylabel('Posterior mean, y_t/n_t','FontSize',myfontsize,'FontName','Times')
%title('Stopping boundary','FontSize',myfontsize,'FontName','Times')
%title('Upper stopping boundary','FontSize',myfontsize,'FontName','Times')
mytitle = strcat('BoundaryYTGreedy','.eps');set(gca,'FontSize',mysmallfontsize);
tmp=axis;
%tmp(1) = max(1,tmp(1));
tmp(1) = 0.9*min(1/scale.gamma./sbvec);
tmp(2) = 1.1*max(1/scale.gamma./sbvec);
%tmp(3)= 10;
%tmp(4) = 1.02*max(up/scale.beta);
axis(tmp);
print('-deps',mytitle);	

