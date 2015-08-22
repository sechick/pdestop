%Test recursion in (w,s) space for sup_S E[ max(W_S/S,0) - 1/S | w_0, s_0]
%problem (the
%normalized problem of maximizing value of decision (pick best of 0 or
%posterior reward) plus sampling cost.  In (y,t) space, we have:
%   sup_T[ max(Y_T/T,0) - c T | y_0/t_0 = \mu_0, t_0] + c t_0,
%for nonanticipating stopping rules T>t_0, and where t = t_0 + n, n is
%number of samples, and y_t = y_0 + \sum_{i=1}^n x_i, X_i ~
%Normal(A,\sigma^2), and prior on A ~ Normal(\mu_0, \sigma_0^2), and
%further t_0 = \sigma^2 / \sigma_0^2.
% updated: 30 Jan 2008

%%% WORK WITH CONSTRUCTOR METHOD

% Set up generic functions for zero discounting
generictermreward=@(wvec,s,p1,p2) max(wvec,0);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
CFApproxValuefunc=@(wvec,s,p1,p2) PDECFKGs(wvec,s);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
upperNoDisc=@(s,p1,p2) CFApproxBoundW(s);
CFfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', CFApproxValuefunc, 'approxmethod', upperNoDisc};
CFscalevec = {'c', 1, 'sigma', 10e5, 'discrate', 0, 'P', 1};
CFparamvec = { 't0', .25, 'tEND', 20000, 'precfactor', 4, 'BaseFileName', 'CF' };
baseparams = { 'online', 0, 'retire', 0, 'DoPlot', 1 };
%figdir Figure\, matdir Matfiles\ UnkVariance 0

% Set up generic functions for positive discounting
CGApproxValuefunc=@(wvec,s,p1,p2) PDECGKGs(wvec,s);   % this is valid terminal reward for discounted rewards, valued in time s currency
upperDisc=@(s,p1,p2) CGApproxBoundW(s);
CGfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', CGApproxValuefunc, 'approxmethod', upperDisc};
CGscalevec = {'c', 0, 'sigma', 10e5, 'discrate', 0.0002, 'P', 1, 'BaseFileName', 'CF' };
CGparamvec = { 't0', .5, 'tEND', 5000, 'precfactor', 4 };

% generic functions when the 'guesses' are still being made regarding the
% upper boundary's approximate value.
Guessdw=0.06; GuessNumW=1500; % these specific values are appropriate for case of P=1, c=1, discrate = 0.0002, for example
upperguessNoClue = [Guessdw GuessNumW]; % guesses for initial dw size, and for number of grid points above and below 0
Guessfunctionset = {'termrewardfunc', generictermreward, 'approxvaluefunc', generictermreward, 'approxmethod', upperguessNoClue};

% Create functions for use with zero discounting
scalevec = CFscalevec; 
paramvec = [CFparamvec, CFfunctionset, baseparams];
[scale, param] = PDEInputConstructor( scalevec, paramvec )
[scale, param] = PDEInputValidator( scale, param )

% Do the computations
tic
[rval, MAXFiles] = PDESolnCompute(scale, param);
toc

% Load in the data structures form those computations
BaseFileName = strcat(param.matdir,param.BaseFileName); % note, we wish to allow loading files by name without having the full solution or the full param: just the name and range of blocks to load
flo = 1;
fhi = MAXFiles;
[rval, cfSoln] = PDESolnLoad(BaseFileName,flo,fhi);
if cfSoln.Header.PDEparam.DoPlot % do a bunch of diagnostics plots, save the eps files
    UtilPlotDiagnostics(cfSoln);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD STUFF %%%%%%%%%%%%%%%%%%%%%%
%SET UP INITIAL PARAMETERS FOR PROBLEM
sigma=100000;     % variance per single sample in (y,t) space
discrate = 0.0002;% zero discount rate per sample for C&F paper (use >0 for C&G paper)
discrate = 0.00;% zero discount rate per sample for C&F paper (use >0 for C&G paper)
if discrate == 0
    c=1;         % cost per single sample in (y,t) space
    t0=0.25;%0.25; %.25
    tEND=30000;%20000; %5000
else
    c=0;
    t0=0.5;%0.25; %.25
    tEND=30000; %25000;20000; %5000
end
m=0;         % value of retirement option in (y,t) space
mu0=0;

APrioriRegret=(sigma/sqrt(t0)) * PsiNorm( abs(mu0-m) / (sigma/sqrt(t0)));

% Set up two data structures to call PDE computation routines.
% 1) Set up parameters of diffusion problem which can be recaled in some
% way. These can be changed without recomputing the solution
PDEscale.c=c;           % cost per single sample in (y,t) space
PDEscale.sigma=sigma;	% standard deviation of a single sample in (y,t) space
PDEscale.P=1;           % multiplier times the mean reward for the adoption decision
PDEscale.discrate = discrate;% zero discount rate per sample for C&F paper (use >0 for C&G paper)
%PDEscale.discrate = 0.00;% zero discount rate per sample for C&F paper (use >0 for C&G paper)
% the following code is a test to see if these values convert ok to scaled
% problem...
PDEscale = PDEScaleStandardize(PDEscale);    % compute alpha, beta, gamma, kappainv for this problem
alpha = PDEscale.alpha;
beta = PDEscale.beta;
gamma = PDEscale.gamma;
kappainv = PDEscale.kappainv;

% 2) Set up some computation specific parameters
%   a) Some involve the terminal reward function for the stopping, which
%   will in general be a function of the scaled posterior means in the
%   variable(vector) wvec, when stopping at time s (valued in time s
%   currency). Normally the parameters p1 and p2 are ignored, but they
%   might be used to (optionally) pass PDEscale (p1) and PDEparam (p2) so as to
%   compute the values. (will especially be useful for the ApproxValuefunc
%   field).
% ** Terminal reward functions: What is expected value of stopping at time
% s with means in wvec
CFTermRewardfunc=@(wvec,s,p1,p2) max(wvec,0);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
CGTermRewardfunc=@(wvec,s,p1,p2) max(wvec,0);   % this is valid terminal reward for discounted rewards, valued in time s currency
% ** Approximate value to go function: if you don't have one for a new
% terminal reward function you might have, then set this to be the same as
% for the terminal reward function. if you have better, e.g. VOI of a fixed
% length sampling policy, a la LL aka KGbeta, or KG*, then it can be used -
% doing so will help to reduce the bias near time 1/(gamma*tEND) when
% starting up the recursion process
'SEC: need to fix the approx value functions'
CFApproxValuefunc=@(wvec,s,p1,p2) PDECFKGs(wvec,s);   % this is valid terminal reward for undiscounted rewards, valued in time s currency
CGApproxValuefunc=@(wvec,s,p1,p2) PDECGKGs(wvec,s);   % this is valid terminal reward for discounted rewards, valued in time s currency
% ** An approximation to the optimal upper stopping boundary. This is a bit of a
% chicken or the egg problem. Once the recursion has been run, you may have
% an approximation - or you might have some analytical results such as from
% the early Chernoff papers. If you have one, pass the function. This will be used 
% for several purposes, including setting the value of dw correctly, and for saving space.
% If you don't have it, you will need to play with setting
%   - dw: bigger to run faster, smaller for more accuracy
%   - NumDw: number of dw in the grid to keep
upperNoDisc=@(s,p1,p2) CFApproxBoundW(s);
upperDisc=@(s,p1,p2) CGApproxBoundW(s);
Guessdw=0.06; GuessNumW=1500; % these specific values are appropriate for case of P=1, c=1, discrate = 0.0002, for example
upperguessNoClue = [Guessdw GuessNumW]; % guesses for initial dw size, and for number of grid points above and below 0
%
if false  % change this to evaluate to true if you want to try the 'ad hoc' computation which does not relay on
    BaseFileName='WildWest'
    PDEparam.termrewardfunc = @(wvec,s,p1,p2) max(wvec,0); % or some other function which is the expected value of stopping at a given time in (w,s) coordinates
    PDEparam.approxvaluefunc = PDEparam.termrewardfunc ; % can use some other function, should be at least as big as termrewardfunction
    PDEparam.approxmethod = upperguessNoClue;
elseif PDEscale.discrate == 0 % set up parameters for undiscounted case of C&F
    % Set up file name for saving results
    BaseFileName='CF';
    PDEparam.termrewardfunc = CFTermRewardfunc;
    PDEparam.approxvaluefunc = CFApproxValuefunc ;
    PDEparam.approxmethod = upperNoDisc;
else % set up parameters for discounted case of C&G
    % Set up file name for saving results
    BaseFileName='CG';
    'SEC: need to validate discounted reward computation'
    PDEparam.termrewardfunc = CGTermRewardfunc;
    PDEparam.approxvaluefunc = CGApproxValuefunc;
    PDEparam.approxmethod = upperDisc;
end
%   b) Some parameters to drive the computations, such as the range of values of n0 to compute, from t0 to tEND,
%       where unknown mean W ~ Normal(mu0, sigma^2/n0), online says if
%       there is online learning, etc.
PDEparam.online = false;% true if online learning is active, false otherwise
PDEparam.t0 = t0;      % lower range of values of t0 to compute, where unknown mean W ~ Normal(mu0, sigma^2/t0)
PDEparam.tEND = tEND;  % upper range of values of t0 to compute, where unknown mean W ~ Normal(mu0, sigma^2/t0)
PDEparam.retire = m;    % value of retirement option: ignored for discrate=0, taken to be alternative for discrate>0
PDEparam.finiteT = false; % false or negative for infinite horizon problem, or a positive finite time horizon. 
% Even if false, need to set tEND to something positve as a 'boundary' condition/time for stopping the PDE solution
% when finitT is true, then the approxvaluefunc is ignored and termrewardfunc is used instead

%   c) Some parameters which guide the computation of the solution and
%   diagnostics during the computation.
if PDEscale.discrate > 0    % precision factor can be smaller for discounted rewards, as larger dw can be used
    PDEparam.precfactor = 10.0;	% increase to increase the precision of the 'dw' grid in normed coordinates. must be at least 1.0.
else
    PDEparam.precfactor = 30.0;	% increase to increase the precision of the 'dw' grid in normed coordinates. must be at least 1.0.
end
PDEparam.DoPlot = true;       % true to turn on diagnostic plots

PDEscale
PDEparam

% Do the computations
tic
[rval, MAXFiles] = PDESolnCompute(BaseFileName, PDEscale, PDEparam);
toc

% Load in the data structures form those computations
flo = 1;
fhi = MAXFiles;
[rval, cfSoln] = PDESolnLoad(BaseFileName,flo,fhi);
if cfSoln.Header.PDEparam.DoPlot % do a bunch of diagnostics plots, save the eps files
    UtilPlotDiagnostics(cfSoln);
end

APrioriRegret=(sigma/sqrt(t0)) * PsiNorm( abs(mu0-m) / (sigma/sqrt(t0)))
tmppp=PDEGetVals(cfSoln,(mu0-m)*beta,1/(t0*gamma))/beta
APrioriRegret-tmppp


% try to check smoothness of value function in these files from one 'block'
% to the next, as a sanity check for numerical stability and the impact of
% the 'ripple' effect
s0 = 1/(gamma * t0)
w0 = mu0*beta
for j=1:(MAXFiles-1)
    j
    svaltotest=cfSoln.Header.lasts(j)
    vala = interp2(cfSoln.Data(j).svec,cfSoln.Data(j).wvec,cfSoln.Data(j).Bwsmatrix,svaltotest,cfSoln.Data(j).wvec)/beta; 
    valb = interp2(cfSoln.Data(j+1).svec,cfSoln.Data(j+1).wvec,cfSoln.Data(j+1).Bwsmatrix,svaltotest,cfSoln.Data(j).wvec)/beta; 
    figure(1111+j);
    plot(cfSoln.Data(j).wvec,vala,'-r',cfSoln.Data(j).wvec,valb,'-.k')
    legend('block j','block j+1');
    figure(1111+MAXFiles+j);
    plot(cfSoln.Data(j).wvec,(vala-valb)./vala,'-',cfSoln.Data(j).upvec(end),0,'x',cfSoln.Data(j).downvec(end),0,'x')
    maxrelerr=max(abs(vala-valb)./vala)
    maxabserr=max(abs(vala-valb))
    legend('relative error: block j - block j+1');
 %   pause;
end
    

% start to reconstruct the values for plotting....
sbvec = cfSoln.Computed.accumsvec;
up = cfSoln.Computed.accumuppper;
lo = cfSoln.Computed.accumlower;
findex = cfSoln.Computed.fileindx;

if PDEparam.DoPlot
    if isa(cfSoln.Header.PDEparam.approxmethod,  'function_handle')
    %    dt = 1/gamma/s0 - 1/gamma/(s0+ds)
    %    dmu=dw/beta
        mysmallfontsize = 14;
        myfontsize = 16;
        figure(17)
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

        figure(18)
        loglog(1/gamma./sbvec,up/beta,'--',1/gamma./sbvec,tstcurv/beta,'-.');set(gca,'FontSize',mysmallfontsize);
        legend('From PDE','Quick Approx.','Location','NorthEast')
        xlabel('Effective number of replications, n_t','FontSize',myfontsize,'FontName','Times'); ylabel('Posterior mean, y_t/n_t','FontSize',myfontsize,'FontName','Times')
        %title('Upper stopping boundary','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat('BoundaryYT','.eps');
        print('-deps',mytitle);	
    end
end


    mwig = PDEparam.retire*beta;
    
    % GET INDEX INTO RIGHT STUCTURE BASED ON TIME
    [a, b] = min(s0 > sbvec);
    ind=findex(b)
    wvec = cfSoln.Data(ind).wvec;
    svec = cfSoln.Data(ind).svec;
    Bwsmatrix = cfSoln.Data(ind).Bwsmatrix;
    Bmu0t0=mwig/beta+interp2(svec,wvec,Bwsmatrix,s0,w0)/beta
    Bw0s0=mwig + interp2(svec,wvec,Bwsmatrix,s0,w0)
    % This is the VOI info - need to add in the stopping to get the full
    % value function
end

PDEGetVals(cfSoln,wvec,sval) % try to get solution for values of w in wvec at time sval, given pde solution cfSoln



%%%%%%%%% convert to standardized version in (w,s) with w brownian motion
%%%%%%%%% in -s time scale
mwig = PDEscale.beta * PDEparam.retire
w0= PDEscale.beta*(mu0-PDEparam.retire)
s0 = 1 / (PDEscale.gamma * PDEparam.t0)
sEND = 1 / (PDEscale.gamma * PDEparam.tEND)

% set up info for plotting stuff
myfontsize=16
mysmallfontsize=14
points = 144*3 %spacing between labels on contours - made so that only one label appears per line


%[ta tb tc]=fminbnd('Bztauapproxc',0,1/gamma/n0base,optimset('TolFun',beta/sigma,'MaxIter',10^4),0,1/gamma/n0base);       % FIXED VERSION
[ta tb tc]=fminbnd('Bztauapproxc',0,s0,optimset('TolFun',beta/sigma,'MaxIter',10^4),0,s0);       % FIXED VERSION
-tb/beta - c/gamma
sig0*PsiNorm(-mu0/sig0)

%ds
%dt = 1/gamma/s0 - 1/gamma/(s0+ds)
%dw
%dmu=dw/beta

%Bmu0t0=m+interp2(svec,wvec,Bwsmatrix,s0,w0)/beta
%Bw0s0=beta*m + interp2(svec,wvec,Bwsmatrix,s0,w0)

% CHUNK5: Generate some plots
% Convert from (w,s) coordinates to (y/t, t) coordinates (NOT (y,t)
% coordinates!!!!)

% REASSEMBLE THE BOUNDARY FROM THE VARIOUS FILES
for ijk=1:MAXFiles
%   reload data
    mymat = strcat(strcat(BaseFileName,int2str(ijk)),'.mat');
%    save(mymat,'svec','wvec','Bwsmatrix','upvec','up1');
    S=load(mymat);
    if ijk==1
        sbvec = S.svec;
        wvec = S.wvec;
        Bws = S.Bwsmatrix;
        up = S.upvec;
        up1 = S.up1;
    else
        sbvec = [sbvec S.svec];
        up = [up S.upvec];
        up1 = [up1 S.up1];
    end
end
figure(17)
tstcurv=ApproxBoundW(sbvec);
plot(sbvec,up,'--',sbvec,tstcurv,'-.');set(gca,'FontSize',mysmallfontsize);
legend('From PDE','Quick Approx.','Location','SouthEast')
xlabel('Reverse time scale,  s=1/(\gamma t)','FontSize',myfontsize,'FontName','Times'); ylabel('Rescaled mean, w','FontSize',myfontsize,'FontName','Times')
%title('Upper stopping boundary','FontSize',myfontsize,'FontName','Times')
mytitle = strcat('BoundaryWS','.eps');
print('-deps',mytitle);	

betastepvec=[1 5 10 20 50 100 200 500 1000];
numbetas=length(betastepvec);
betamatrix=zeros(numbetas,length(up));

figure(18)
loglog(1/gamma./sbvec,up/beta,'--',1/gamma./sbvec,tstcurv/beta,'-.');set(gca,'FontSize',mysmallfontsize);
legend('From PDE','Quick Approx.','Location','NorthEast')
xlabel('Effective number of replications, n_t','FontSize',myfontsize,'FontName','Times'); ylabel('Posterior mean, y_t/n_t','FontSize',myfontsize,'FontName','Times')
%title('Upper stopping boundary','FontSize',myfontsize,'FontName','Times')
mytitle = strcat('BoundaryYT','.eps');
print('-deps',mytitle);	

figure(19)
tstcurv=ApproxBoundW(sbvec);
plot(up./tstcurv)

% Try to get one-step estimate of boundary
for j=1:numbetas;
    ybndup1step = up/beta; %initialize vector
    tvec = 1/gamma./sbvec;
    ytinit = ybndup1step(length(ybndup1step))*tvec(length(ybndup1step));
    betamatrix(j,:) = ytinit;
    nrepslookahead=betastepvec(j)
    for i=length(ybndup1step):-1:1
        [ytinit, fval, exitflag]=fzero(@PsiNormRepsRoot,ytinit,optimset('TolX',1e-8),sigma,tvec(i),nrepslookahead,nrepslookahead*c);
        ybndup1step(i)=ytinit;
        ytinit=ytinit*0.99;
    end
    betamatrix(j,:)=ybndup1step;
end
kgstar=max(betamatrix);

figure(19)
%loglog(1/gamma./sbvec,up/beta,'--',1/gamma./sbvec,tstcurv/beta,'-.',tvec,ybndup1step./tvec,':o',tvec,ybndupNstep./tvec,'-x')
%legend('From PDE','Quick Approx.','One-step lookahead','N-step lookahead','Location','NorthEast')
loglog(1/gamma./sbvec,up/beta,'-',1/gamma./sbvec,kgstar./tvec,'--',tvec,betamatrix(1,:)./tvec,':',tvec,betamatrix(3,:)./tvec,'.',tvec,betamatrix(6,:)./tvec,'-.')
legend('PDE','KG_*','KG_1 = 1-step lookahead','KG_{10} = 10-step lookahead','KG_{100} =100-step lookahead','Location','SouthWest')
%loglog(1/gamma./sbvec,tstcurv/beta,'-.',tvec,ybndup1step./tvec,':')
%legend('PDE or Quick Approx.','One-step lookahead','Location','NorthEast')
xlabel('Effective number of replications, n_t','FontSize',myfontsize,'FontName','Times'); ylabel('Posterior mean, y_t/n_t','FontSize',myfontsize,'FontName','Times')
%title('Stopping boundary','FontSize',myfontsize,'FontName','Times')
%title('Upper stopping boundary','FontSize',myfontsize,'FontName','Times')
mytitle = strcat('BoundaryYTGreedy','.eps');set(gca,'FontSize',mysmallfontsize);
tmp=axis;
tmp(1) = max(1,tmp(1));
tmp(2) = 10^4; max(1/gamma./sbvec);
tmp(3)= 10;
tmp(4) = 5*10^5; 1.02*max(up/beta);
axis(tmp);
print('-deps',mytitle);	



% Try to get one-step estimate of boundary
ybndup1step = up/beta; %initialize vector
tvec = 1/gamma./sbvec;
ytinit = ybndup1step(length(ybndup1step))*tvec(length(ybndup1step));
nrepslookahead=1;
for i=length(ybndup1step):-1:1
    [ytinit, fval, exitflag]=fzero(@PsiNormRepsRoot,ytinit,optimset('TolX',1e-8),sigma,tvec(i),nrepslookahead,nrepslookahead*c);
    ybndup1step(i)=ytinit;
    ytinit=ytinit*0.99;
end

% Try to get N-step estimate of boundary
ybndupNstep = ybndup1step; %initialize vector
ytinit = ybndupNstep(length(ybndupNstep))*tvec(length(ybndupNstep)) * 1.01;
nrepslookahead=10;
for i=length(ybndupNstep):-1:1
    [ytinit, fval, exitflag]=fzero(@PsiNormRepsRoot,ytinit,optimset('TolX',1e-8),sigma,tvec(i),nrepslookahead,nrepslookahead*c);
    ybndupNstep(i)=ytinit;
    ytinit=ytinit*1.01;
end

% Try to get another N-step estimate of boundary
ybndupNbisstep = ybndup1step; %initialize vector
ytinit = ybndupNstep(length(ybndupNstep))*tvec(length(ybndupNstep)) * 1.01;
nrepslookahead=100;
for i=length(ybndupNstep):-1:1
    [ytinit, fval, exitflag]=fzero(@PsiNormRepsRoot,ytinit,optimset('TolX',1e-8),sigma,tvec(i),nrepslookahead,nrepslookahead*c);
    ybndupNbisstep(i)=ytinit;
    ytinit=ytinit*1.01;
end


tmpvec=ybndupNstep-ybndup1step;
max(tmpvec)
min(tmpvec)

figure(19)
%loglog(1/gamma./sbvec,up/beta,'--',1/gamma./sbvec,tstcurv/beta,'-.',tvec,ybndup1step./tvec,':o',tvec,ybndupNstep./tvec,'-x')
%legend('From PDE','Quick Approx.','One-step lookahead','N-step lookahead','Location','NorthEast')
loglog(1/gamma./sbvec,up/beta,'-',tvec,ybndup1step./tvec,':',tvec,ybndupNstep./tvec,'-.',tvec,ybndupNbisstep./tvec,'--')
legend('PDE','\beta=  1-step lookahead','\beta= 10-step lookahead','\beta=100-step lookahead','Location','NorthEast')
%loglog(1/gamma./sbvec,tstcurv/beta,'-.',tvec,ybndup1step./tvec,':')
%legend('PDE or Quick Approx.','One-step lookahead','Location','NorthEast')
xlabel('Effective number of replications, n_t','FontSize',myfontsize,'FontName','Times'); ylabel('Posterior mean, y_t/n_t','FontSize',myfontsize,'FontName','Times')
%title('Stopping boundary','FontSize',myfontsize,'FontName','Times')
%title('Upper stopping boundary','FontSize',myfontsize,'FontName','Times')
mytitle = strcat('BoundaryYTGreedy','.eps');
print('-deps',mytitle);	

