function [rval, numfiles] = PDESolnCompute(fName, PDEscalein, PDEparam)
% CFSolvePDE: Solves free boundary problem in scaled coordinates for the
% sequential sampling for simulation optimization. 
%
% For sample calling conventions, including how to set up the fields of the
% structures PDEscale and PDEparam, please see TestSolvePDE.m.
%
% INPUTS: 
%   fName: base name for files to hold the data structures with the
%       solutions of the free boundary PDE.
%   PDEscale: contains basic parameters of the problem which might be
%      changed without need to recompute the PDE solution (c, sigma, 
%      discrate, P). The solution computation will compute the relevant 
%      standardized values (alpha, beta, gamma, kappainv)
%   PDEparam: contains parameters of the problem which might entail
%      recomputation of the solution if they are changed (online learning), 
%      or which specify details of how the computation is to be done (e.g.
%      terminal reward function, approximation to the value function, upper
%      and lower ranges for t in (y/t, t) space for specifying probability
%      distrubitions, granularity of matrix computations, flag for 
%      diagnostic plots, etc.
%       - PDEparam.precfactor = 4.0;	% increase to increase the precision of the 'dw' grid in normed coordinates. must be at least 1.0.
%       - PDEparam.DoPlot = true;       % true to turn on diagnostic plots
%
% OUTPUTS:
%   rval: true if completed ok, false if there were errors
%   numfiles: number of files output
%   Side effect: files of name fName<n>.mat for n=1,2,...,numfiles
%       which contain solutions for PDE in standardized coordinates 
%       Contents include: upper and lower boundaries, standardized: reward
%       function, PCS, expected number of samples to stop time, ...
%       also creates BASEFILName0.mat, with summary information for files.
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
%

% standardize the scaling again
PDEscale = PDEScaleStandardize(PDEscalein);

% copy some parameter values to local
%alpha = PDEscale.alpha;
beta = PDEscale.beta;
gamma = PDEscale.gamma;
%kappainv = PDEscale.kappainv;
termreward = PDEparam.termrewardfunc;   % this function should return the reward for stopping at the specified time
if PDEparam.finiteT
    termreward = PDEparam.termrewardfunc;   % this function should return the reward for stopping at the specified time
else
    rewardfunc = PDEparam.approxvaluefunc;  % this function should return the best approximation for the reward to go, or can be termreward if the original problem is a finite time horizon problem
end
approxmeth = PDEparam.approxmethod;  % this function should return the best approximation for the reward to go, or can be termreward if the original problem is a finite time horizon problem

s0 = 1 / (gamma * PDEparam.t0);
sEND = 1 / (gamma * PDEparam.tEND);

mwig = beta * PDEparam.retire;
%muwig0 = beta*PDEparam.mu0;

%%% Do some input validation
rval = 1;
numfiles = 0;
if length(fName) < 1
    fName = 'PDEData';
end
if (min(s0,sEND) <= 0) || (s0 <= sEND)
    warning('must have s0 > sEND > 0');
    rval = 0;
end
isDisc = (PDEscale.discrate > 0);
isOffline = ~PDEparam.online;

if isDisc
    warning('computation for discounted rewards still being validated');
%    rval = 0;
end
if ~isOffline
    warning('computation for online rewards not yet implemented');
    rval = 0;
end
if PDEparam.precfactor < 1.0
    warning('Value of precfactor should be at least 1.0.  Changing precfactor to 1.0');
	PDEparam.precfactor = 1.0;
end
if ~rval
    return;
end

myfontsize=16;
mysmallfontsize=14;
points = 144*3; %spacing between labels on contours - made so that only one label appears per line

% This code follows the idea of Chernoff & Petkau, and of Brezzi and Lai,
% papers, except that we use a trinomial tree rather than binomial tree.

% Start setting up computational variables
if PDEparam.finiteT
    sinit = sEND;      % start the calculation of the free boundary at values beyond sEND (smaller than sEND) in reverse time scale
else
    sinit = 0.95*sEND;      % start the calculation of the free boundary at values beyond sEND (smaller than sEND) in reverse time scale
end
if isa(approxmeth, 'function_handle')
%    if ~isDisc % for discounted case, want to have precfactor>1 grid points above 0 where expected reward for continuing is at least 0
%        dw = sinit^2/(20/9)/PDEparam.precfactor;  % need to divide be something bigger than 2 in order to get a valid recursion when C&F model is used
%    else       % for discounted case, want at least precfactor>1 grid points below the asymptotic approx for boundary for small sinit
%        dw = (min(0.2,sinit)/sqrt(2))/PDEparam.precfactor/(10/9);      % need to divide be something bigger than 2 in order to get a valid recursion when C&F model is used
%    end
    dw = approxmeth(sinit,PDEscale,PDEparam) / PDEparam.precfactor / (10/9);    % the extra 10/9 is a fudge factor for hoping for numerical error improvement
else
    dw = approxmeth(1);
end
ds = dw^2 * 2 / 3;          % in trinomial tree, equal probs of going up, straight or down implies this equation to get correct variance of reverse time brownian motion
SAVEEVERY=160;        % every how many iterations do we keep the PDE values, for storage in matrix and files?
NUMSAVESPERITER=100;  % number of times we save the PDE values, for storage in matrix and files?
Numds=SAVEEVERY*NUMSAVESPERITER;   % number of time steps per grid, % 8000 or 16000 for example, 
FRACTOKEEP = 0.9; % use this to remember the state at some earlier time: 0.9 or 085 should get of ripples in many cases, make sure it is smaller than
if (FRACTOKEEP > (1-2/SAVEEVERY)) | (FRACTOKEEP < 0.5)
    warning('keep FRACTOKEEP between 0.5 and 1 - 2/SAVEEVERY, please insure SAVEEVERY big enough');
end
myeps = min(1,beta)*min(10^-12, ds/1000);    % a small value, for use in checking if things are close to 0 or not

if PDEparam.DoPlot
    s0
    sEND
end


%%%%%%%%%%%% START CLIPPING INTO MATLAB ALL OF THIS, CHUNK BY CHUNK
% CHUNK1: Set up the Grid

% now figure out what is the right number of dw values which are needed
% along the way - knowing that dw and ds will be scaled up, and the
% boundaries are roughly approximated using the asymptotic
% approximations...
THRESHSCALE = 100;          % set to be a bit bigger than the threshold value for 's' in the scaled approx for optimal boundary
numchunks = 1 + ceil( log2( 1.1*(min(THRESHSCALE,s0)-sinit)/(Numds*ds) ) ); % guess as to how many iterations are needed to generate all the files
spowvec = 4.^(0:(numchunks-1)); % compute some vectors to say roughly how many times the initial 'ds' are the s-scale steps in the i=1,...,numchunks iterations
wpowvec = 2.^(0:(numchunks-1)); % compute some vectors to say roughly how many times the initial 'dw' are the w-scale steps in the i=1,...,numchunks iterations
smaxvec = sinit+Numds*ds*cumsum(spowvec); % find the 
if isa(approxmeth, 'function_handle')
    wmaxvec=approxmeth(smaxvec,PDEscale,PDEparam);      % approximate the maximum value of w for which the grid must account
    dwvec = dw*wpowvec;
    %dsvec = ds*spowvec;
    ceilfactor=2.0;
    if isDisc   % get an approximation for the range of the upper and lower boundaries
        lowceilfactor = 2.0*ceilfactor;   % make fudge factor bigger for discounted case, as the range for a lower boundary is wider...
    end
    ratiovec = wmaxvec ./ dwvec;    % figure roughly how many dw from 0 to upper boundary...
    bigw = ceil(max(ceilfactor*ratiovec)); % hit that figure with a multiplier in order to make sure there is a margin of error, especially below boundary
else  % if we don't have a function to compute, approximately, the upper stopping boundary, then take it from user specified input
    bigw = max(2*PDEparam.precfactor,ceil(approxmeth(2)));
end

% Now, can initialize the vector of w values for the recursion, and compute
% how big it is.
wvec = dw*(-bigw:bigw);
[wvsize] = length(wvec);

% Set up the initial conditions for the recursions
scur = sinit
% vectors for value function
Cinitvec = rewardfunc(wvec,scur,PDEscale, PDEparam); % Approximation of value to go, given stopping at the discretized time sinit.
Cin=Cinitvec;    % initialize terminal condition assuming immediate stopping.
Cout=Cin;           % only to initialize
% vectors for EN = expected number of samples
ENCin = 0*Cinitvec; % a priori, assume 0 more samples in expectation at terminal time
ENCout = ENCin; % a priori, assume 0 more samples in expectation at terminal time
% structures for PCS
EPCSCout = 0*Cinitvec; % a priori, assume 0 more samples in expectation at terminal time
EPCSCin = 0*Cinitvec; % a priori, assume 0 more samples in expectation at terminal time

% CHUNK2: Defensive coding section, to check if ds and sinit are ok in size relative to each other.
% Try to iterate to find a value of s such that the continuation set
% at least includes the value w_s=0.  Note that the discretization may make
% it so taht the approximation puts (w_s,s) in the stopping set even though
% it is in theory in the continuation set... so if it is hard to find such
% an s then ds is likely to be too large.  
% If, in your application this runs for a lot of iterations, then probably
% it is worth deriving the appropriate size for dw and ds, relative to the
% initial starting position, sinit, as small sinit might impose some
% constraints. 
% For the CF code, the derivations should make this so it only needs to
% loop the minimum number of times.
middlecin1 = Cin(bigw-3:bigw+3)-Cinitvec(bigw-3:bigw+3)
minindx=1; counter = 1;
while minindx==1
    % do the iteration
    sout = scur+ds;
    Cinitvec = termreward(wvec,sout,PDEscale, PDEparam);  % compute reward from stopping
    if ~isDisc % for undiscounted rewards, do the update, noting that the top and bottom can be 'tied' to be 'stopped' 
        Cout(2:(wvsize-1)) = ( Cin(1:(wvsize-2)) + Cin(2:(wvsize-1)) + Cin(3:wvsize)) / 3  - (1/scur-1/sout);
        Cout(1) = Cinitvec(1);              % following conditions in sure stopping at top and at bottom of vector
        Cout(wvsize) = Cinitvec(wvsize);
    else  % with discounting, it might be possible to have no lower stopping region, therefore we need to handle the update for Cout(1) in special way
        % we use a surrogate for the reward of doing dw below the stopping boundary... 
        Cout(2:(wvsize-1)) =  exp(- (1/scur-1/sout)) * ( Cin(1:(wvsize-2)) + Cin(2:(wvsize-1)) + Cin(3:wvsize)) / 3 ;
        Cout(1) = exp(- (1/scur-1/sout)) * (termreward(min(wvec)-dw,sout,PDEscale, PDEparam) + Cin(1) + Cin(2) ) / 3 ;   % following conditions in sure stopping at top and at bottom of vector
        Cout(wvsize) = Cinitvec(wvsize);
    end
    Cout=max(Cinitvec,Cout);
    Cin=Cout;                           % iterate to obtain new initial condition

    scur=sout;
    
    stopped = ( abs(Cout - Cinitvec) < myeps ); % check if we have stopped or not.
    [~, minindx] = min(stopped);             % get index to first grid point above lower stopping boundary
    [~, maxindx] = min(flipdim(stopped,2));  % get index to first grid point below upper stopping boundary
    maxindx = wvsize + 1 - maxindx;
    if (counter < SAVEEVERY) & (~PDEparam.finiteT)  % force a few iterations of this recursion in order to 'prime' the initial conditions and to avoid 'edge' effects from the initializaiton of the terminal condition
%    if (counter < 20*PDEparam.precfactor)  % force a few iterations of this recursion in order to 'prime' the initial conditions and to avoid 'edge' effects from the initializaiton of the terminal condition
        minindx = 1;
    end
    counter = counter + 1;
end 
if PDEparam.DoPlot
    dw
    bigw
    scur
%    numinitchecks = (scur - sinit)/ds
%    scurminindxmaxindx = [scur minindx maxindx]
%    middlecin = Cin(minindx-3:maxindx+3)-Cinitvec(minindx-3:maxindx+3)
    middlecin = Cin(maxindx-3-PDEparam.precfactor:maxindx+3)-Cinitvec(maxindx-3-PDEparam.precfactor:maxindx+3)
end

% CHUNK3: Now that we have a point (w=0,s) that is in the continuation set, the
% next set of code initializes vectors/matrices to start collecting data
% about the continuation set and the value B(.,.) of the value function
% set up data structures:
svec=zeros(1,1+NUMSAVESPERITER); %1+floor((Numds-1)/SAVEEVERY));
%Vwsmatrix=zeros(wvsize,1+floor((Numds-1)/SAVEEVERY));
Bwsmatrix=zeros(wvsize,1+NUMSAVESPERITER); %floor((Numds-1)/SAVEEVERY));
ENwsmatrix=0*Bwsmatrix;
EPCSwsmatrix=0*Bwsmatrix;
up1=zeros(size(svec));
upvec=zeros(size(svec));
down1=zeros(size(svec));
downvec=zeros(size(svec));

% CHUNK4: do the iteration for the value function, only saving the data structure
% every so often, to reduce data storage requirements.
ijk=0;
while (sout < s0) %&& (wmax ~= wvec(maxindx))        % iterate until the largest value of s0 that is needed is covered.
    ijk=ijk+1;          % counter for iterations
    firsts(ijk) = scur+ds;     % remember smallest value of s for which this block computes statistics
    for i=1:Numds+1 % run one extra, so first and last iterations of data are saved (see the mod(...) == 1 below)
        % do the iteration of the first order PDE solution
        sout = scur+ds;
        Cinitvec = termreward(wvec,sout,PDEscale, PDEparam);  % compute reward for stopping at time sout
        
        if ~isDisc % for undiscounted rewards, do the update, noting that the top and bottom can be 'tied' to be 'stopped' 
            Cout(2:(wvsize-1)) = ( Cin(1:(wvsize-2)) + Cin(2:(wvsize-1)) + Cin(3:wvsize)) / 3  - (1/scur-1/sout);
            Cout(1) = Cinitvec(1);              % following conditions in sure stopping at top and at bottom of vector
            Cout(wvsize) = Cinitvec(wvsize);
        else  % with discounting, it might be possible to have no lower stopping region, therefore we need to handle the update for Cout(1) in special way
            % we use a surrogate for the reward of doing dw below the stopping boundary... 
            Cout(2:(wvsize-1)) =  exp(- (1/scur-1/sout)) * ( Cin(1:(wvsize-2)) + Cin(2:(wvsize-1)) + Cin(3:wvsize)) / 3 ;
            Cout(1) = exp(- (1/scur-1/sout)) * (rewardfunc(min(wvec)-dw,sout,PDEscale, PDEparam) + Cin(1) + Cin(2) ) / 3 ;              % following conditions in sure stopping at top and at bottom of vector
            Cout(wvsize) = Cinitvec(wvsize);
        end
        Cout=max(Cinitvec,Cout);
        Cin=Cout;                           % iterate to obtain new initial condition

        % try to compute the stopping boundaries
        stopped = ( abs(Cout - Cinitvec) < myeps );
        [~, minindx] = min(stopped);             % get index to first grid point above lower stopping boundary
        [~, maxindx] = min(flipdim(stopped,2));  % get index to first grid point below upper stopping boundary
        maxindx = wvsize + 1 - maxindx;         %    [scur minindx maxindx];

        % FIX: update expected number of samples til stopping time: to be
        % debugged
        ENCout(2:(wvsize-1)) = ( ENCin(1:(wvsize-2)) + ENCin(2:(wvsize-1)) + ENCin(3:wvsize)) / 3  + (1/scur-1/sout);
        ENCout(stopped) = 0;
        ENCin = ENCout; 

        % FIX: compute PCS: to be debuged. 
        EPCSCout(2:(wvsize-1)) = ( EPCSCin(1:(wvsize-2)) + EPCSCin(2:(wvsize-1)) + EPCSCin(3:wvsize)) / 3;
        EPCSCout(stopped) = normcdf(abs(wvec(stopped))./sqrt(sout),0,1); % put z-statistic in stopping zone, which is prob that a normal with mean w and std s is > 0
        EPCSCin = EPCSCout; 

        % update time
        scur=sout;

        % check to see if we should remember the value of the boundary and
        % value function
        if mod(i,SAVEEVERY)==1   % by not saving at first iteration, 
            ind=1+floor(i/SAVEEVERY);
            svec(ind) = scur;                   % store value of s values
    %        Vwsmatrix(:,ind)=Cin;               % store value function approximation
            Bwsmatrix(:,ind)=Cin-Cinitvec;      % store value of information above stopping immediately
            ENwsmatrix(:,ind)=ENCin;            % store value function approximation, above max(w_s,s)
            EPCSwsmatrix(:,ind)=EPCSCin;             % store value function approximation, above max(w_s,s)
            D0 = Cin(maxindx) - Cinitvec(maxindx);          % Examine first two steps not in boundary (Chernoff and Petkau)
            D1 = Cin(maxindx-1) - Cinitvec(maxindx-1);
            up1(ind)=wvec(maxindx);             % the 'raw' value for the boundary of the continuation set
            upvec(ind) = up1(ind) + sqrt(ds) * abs(D1 / (2*D1 - 4*D0)); % the 'corrected' value for the boundary
            D0 = Cin(minindx) - Cinitvec(minindx);          % Examine first two steps not in boundary (Chernoff and Petkau)
            D1 = Cin(minindx+1) - Cinitvec(minindx+1);
            down1(ind)=wvec(minindx);
            downvec(ind) = down1(ind) - sqrt(ds) * abs(D1 / (2*D1 - 4*D0));
    %        % Here we don't compute down boundary to save computational time, since by symmetry
    %        properties it is the negative of the upper boundary
            if (i-2) <= floor(FRACTOKEEP*Numds) 
                lastsaveds=scur;
                lastsavedCin=Cin;
                lastsavedENCin=ENCin;
                lastsavedEPCSCin=EPCSCin;
            end
        end
    end
    lasts(ijk) = scur;     % remember largest value of s for which this block computes statistics
    lasthelds(ijk) = lastsaveds; % remember the value of s for restarting the next block
    
    % save results to file
    mymat = strcat(fName,int2str(ijk),'.mat');
% save(mymat,'svec','wvec','Vwsmatrix','Bwsmatrix','ENwsmatrix','EPCSwsmatrix','upvec','up1','downvec','down1');
% should be able to reconstruct Vwsmatrix from Bwsmatrix and PDEscale, PDEparam.
    save(mymat,'Bwsmatrix','ENwsmatrix','EPCSwsmatrix','svec','wvec','upvec','downvec');

    % if diagnostic plots are desired, plot them
    if PDEparam.DoPlot
        figdir = 'Figure\';
        if ~isdir(figdir) 
            mkdir(figdir);
        end
        
        myfontsize=16;
        mysmallfontsize=14;
        points = 144*3; %spacing between labels on contours - made so that only one label appears per line

        ijk
        maxup1upwbiaswvec = [max(up1) max(upvec) max(wvec)]
        % First, plot a contour plot of the benefit (above 0) of continuing, plus a
        % boundary of the upper and lower continuation set
        figure(20+ijk)
        hold off
        [C, h]=contour(svec,wvec,Bwsmatrix); clabel(C,h,'FontSize',mysmallfontsize,'FontName','Times','LabelSpacing',points);
        hold on
        plot(svec,upvec,'--',svec,downvec,'--');set(gca,'FontSize',mysmallfontsize);
        if ijk>1
            plot([ lasts(ijk-1) lasts(ijk-1) ] ,[min(wvec) max(wvec)],'-.r');
        end
        plot([ firsts(ijk) firsts(ijk) ] ,[min(wvec) max(wvec)],'-.g');
        plot([ lasthelds(ijk) lasthelds(ijk) ] ,[min(wvec) max(wvec)],'-.k');
        tmp=axis;tmp(4)=max(10*dw,1.2*max(max(up1),max(upvec)));tmp(3)=1.2*min(min(down1),min(downvec));axis(tmp);
        xlabel('Reverse time scale, s','FontSize',myfontsize,'FontName','Times'); ylabel('Scaled mean, w_s','FontSize',myfontsize,'FontName','Times')
        title('Stdized E[value of continuing over stopping | (w_s, s)]','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigContWS',int2str(ijk),'.eps');
        print('-deps',mytitle);	

        figure(60+ijk)
        hold off
        [C, h]=contour(1/gamma./svec,wvec/beta,Bwsmatrix/beta); clabel(C,h,'FontSize',mysmallfontsize,'FontName','Times','LabelSpacing',points);
        hold on
        plot(1/gamma./svec,upvec/beta,'--',1/gamma./svec,downvec/beta,'--');set(gca,'FontSize',mysmallfontsize);
        tmp=axis;tmp(4)=1.2*max(max(10*dw,max(up1))/beta,max(upvec)/beta);tmp(3)=1.2*min(min(down1),min(downvec))/beta;axis(tmp);
        xlabel('Effective number of samples, n_t','FontSize',myfontsize,'FontName','Times'); ylabel('Posterior mean, y_t/n_t','FontSize',myfontsize,'FontName','Times')
        title('E[value of continuing over stopping | (y_t/n_t, n_t)]','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigContYT',int2str(ijk),'.eps');
        print('-deps',mytitle);	

        figure(100+ijk)
        hold off
        plot(svec,up1,'-.',svec,upvec,'--');set(gca,'FontSize',mysmallfontsize);
        hold on
        if isa(approxmeth, 'function_handle')
            theoryvec=approxmeth(svec,PDEscale,PDEparam);   
            plot(svec,theoryvec,'-');
            legend('not adjusted','bias adjusted','theory bound');
        else
            legend('not adjusted','bias adjusted');
        end
        plot(svec,down1,'-.',svec,downvec,'--');set(gca,'FontSize',mysmallfontsize);
        tmp=axis;tmp(4)=1.2*max(max(10*dw,max(up1)),max(upvec));tmp(3)=1.2*min(min(down1),min(downvec));axis(tmp);
        xlabel('Reverse time scale, s','FontSize',myfontsize,'FontName','Times'); ylabel('Scaled mean, w_s','FontSize',myfontsize,'FontName','Times')
        title('Stopping boundaries in (w_s,s) scale','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigContCPbias',int2str(ijk),'.eps');
        print('-deps',mytitle);	

        %        figure(140+ijk)
        figure(140)
        hold off
        plot(svec,-downvec,'-.',svec,upvec,'--');set(gca,'FontSize',mysmallfontsize);
        hold on
        if isa(approxmeth, 'function_handle')
            theoryvec=approxmeth(svec,PDEscale,PDEparam);   
            plot(svec,theoryvec,'-');
            legend('- downvec','upvec','theory bound');
        else
            legend('- downvec','upvec');
        end
        tmp=axis;tmp(4)=1.2*max(max(10*dw,max(-downvec)),max(upvec));tmp(3)=0;axis(tmp);
        xlabel('Reverse time scale, s','FontSize',myfontsize,'FontName','Times'); ylabel('Scaled mean, w_s','FontSize',myfontsize,'FontName','Times')
        title('Compare upper bound with -lower bound, (w_s, s) coord','FontSize',myfontsize,'FontName','Times')
        %        mytitle = strcat(figdir,fName,'FigDiffUpDownWS',int2str(ijk),'.eps');
        %        print('-deps',mytitle);	

        %        figure(180+ijk)
        figure(180)
        hold off
        plot(wvec(:),Cinitvec(:),'-.',wvec(:),Cin(:),'--');set(gca,'FontSize',mysmallfontsize);
        hold on; 
        plot(downvec(ind),0,'x',upvec(ind),0,'x');set(gca,'FontSize',mysmallfontsize);
        % tmp=axis;tmp(4)=1.2*max(max(10*dw,max(up1))/beta,max(upvec)/beta);tmp(3)=0;axis(tmp);
        xlabel('wvec','FontSize',myfontsize,'FontName','Times'); ylabel('V(s,w)','FontSize',myfontsize,'FontName','Times')
        title('Value function in (w,s) coords','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigVWS',int2str(ijk),'.eps');
        print('-deps',mytitle);	

        %        figure(220+ijk)
        figure(220)
        hold off
        plot(wvec(:)/beta,ENCin(:)/gamma,'-');set(gca,'FontSize',mysmallfontsize);
        hold on; 
        plot(downvec(ind)/beta,0,'x',upvec(ind)/beta,0,'x');set(gca,'FontSize',mysmallfontsize);
        tmp=axis;tmp(4)=1.2*max(ENCin)/gamma;tmp(3)=-1;axis(tmp);
        xlabel('Posterior mean y_t/n_t','FontSize',myfontsize,'FontName','Times'); ylabel('E[num samples]','FontSize',myfontsize,'FontName','Times')
        title('E[num samples] in (y_t/n_t,n_t) coords','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigENYT',int2str(ijk),'.eps');
        print('-deps',mytitle);	

        %        figure(260+ijk)
        figure(260)
        hold off
        plot(wvec(:),1-EPCSCin(:),'-',wvec(:),1 - normcdf(abs(wvec(:))/sqrt(scur),0,1),'-.');set(gca,'FontSize',mysmallfontsize);
        hold on; 
        plot(downvec(ind),0,'x',upvec(ind),0,'x');set(gca,'FontSize',mysmallfontsize);
        %        tmp=axis;tmp(4)=1.2*(1-max(EPCSCin));tmp(3)=0;axis(tmp);
        xlabel('wvec','FontSize',myfontsize,'FontName','Times'); ylabel('1-E[PCS]','FontSize',myfontsize,'FontName','Times')
        title('1-E[PCS given (w, s)]','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigEPCSWS',int2str(ijk),'.eps');
        print('-deps',mytitle);	

        figure(300)
        hold off
        semilogy(wvec(:),Cinitvec(:),'-.',wvec(:),Cin(:),'--');set(gca,'FontSize',mysmallfontsize);
        hold on; 
        legend('Cinitvec','Cin');
        semilogy(downvec(ind),1,'x',upvec(ind),1,'x');set(gca,'FontSize',mysmallfontsize);
        % tmp=axis;tmp(4)=1.2*max(max(10*dw,max(up1))/beta,max(upvec)/beta);tmp(3)=0;axis(tmp);
        xlabel('wvec','FontSize',myfontsize,'FontName','Times'); ylabel('V(s,w)','FontSize',myfontsize,'FontName','Times')
        title('Value function in (w,s) coords','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigLogVWS',int2str(ijk),'.eps');
        print('-deps',mytitle);	

        WidthContinInw = upvec(ind) - downvec(ind)
    end  % diagnostic plot routines

    % double check boundaries to see if bias corrected bound exceeds max w
    % value in calculation or not. if there is an error here, then need to
    % set bigw to a larger value
    if (max(max(up1),max(upvec)) >= max(wvec))
        error('error: need bigger initial number of grid points for w direction');
        return;
    end

    %   now reset the ds and dw values to scale up an order of magnitude
    ds = 4*ds;      % quadruple the time step
    dw = 2*dw;       % double the space step, but keep the same number of steps
    wvec=2*wvec;    % change the w-space vector, so range is twice as wide with same number of steps
    
    % Can't just restart the recursion - this would cause a 'ripple' due to
    % the change in discretization: so what we do is jump back to an
    % 'earlier' time, with the new discretization, and continu the
    % regression for some previously computed values of s, so that the
    % 'ripple' goes away and there is smoothness in the output again.
    % 1) first, reload
    scur=lastsaveds;
    Cinitvec=termreward(wvec,scur,PDEscale, PDEparam); % compute reward for stopping at time scur
    Cin=Cinitvec;   % sets boundary conditions for values above and below.
    midindx=(wvsize+1)/2;
    cntindx=floor((midindx-1)/2);
    for xyz=1:cntindx
        Cin(midindx+xyz)=lastsavedCin(midindx+2*xyz);
        Cin(midindx-xyz)=lastsavedCin(midindx-2*xyz);
        ENCin(midindx+xyz)=lastsavedENCin(midindx+2*xyz);
        ENCin(midindx-xyz)=lastsavedENCin(midindx-2*xyz);
        EPCSCin(midindx+xyz)=lastsavedEPCSCin(midindx+2*xyz);
        EPCSCin(midindx-xyz)=lastsavedEPCSCin(midindx-2*xyz);
    end
    Cin(midindx)=lastsavedCin(midindx);
	ENCin(midindx)=lastsavedENCin(midindx);
    EPCSCin(midindx)=lastsavedEPCSCin(midindx);
%    plot(wvec,Cin,'--')
%    pause
    hold off
end
numfiles=ijk;

% now save header file: store info related to computation of results.
ijk=0;
StartFileVal = 1;           % lower index of valid file values
EndFileVal = numfiles;      % upper index of valid file values
TimeStamp = clock;
mymat = strcat(fName,int2str(ijk),'.mat');
save(mymat,'fName', 'TimeStamp','StartFileVal','EndFileVal','PDEscale','PDEparam','lasts','firsts','lasthelds','myeps');


end