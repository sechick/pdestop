function [rval, numfiles] = PDESolnCompute(PDEscalein, PDEparam)
% CFSolvePDE: Solves free boundary problem in scaled coordinates for the
% sequential sampling for simulation optimization. 
%
% For sample calling conventions, including how to set up the fields of the
% structures PDEscale and PDEparam, please see TestSolvePDE.m.
%
% INPUTS: 
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
%       also creates fBame0.mat, with summary information for files, where
%       fName is determined by PDEparam field.
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

figdir = PDEparam.figdir;
matdir = PDEparam.matdir;
fName = PDEparam.BaseFileName;
%
%figdir = 'Figure\';
figsave = true;

% create the directory for the figures if it does not exist already and the figures are
% to be saved
if ~isdir(figdir) && figsave
    mkdir(figdir);
end

% copy some parameter values to local
%alpha = PDEscale.alpha;
%beta = PDEscale.beta;
gamma = PDEscale.gamma;
%kappainv = PDEscale.kappainv;
termreward = PDEparam.termrewardfunc;   % this function should return the reward for stopping at the specified time
if PDEparam.finiteT
    rewardfunc = PDEparam.termrewardfunc;   % this function should return the reward for stopping at the specified time
else
    rewardfunc = PDEparam.approxvaluefunc;  % this function should return the best approximation for the reward to go, or can be termreward if the original problem is a finite time horizon problem
end
approxmeth = PDEparam.approxmethod;  % this function should return the best approximation for the reward to go, or can be termreward if the original problem is a finite time horizon problem

s0 = 1 / (gamma * PDEparam.t0);
sEND = 1 / (gamma * PDEparam.tEND);

%mwig = beta * PDEparam.retire;
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

if ~isOffline
    warning('computation for online rewards not yet implemented');
    rval = 0;
end
if PDEparam.precfactor < 2.0
    warning('Value of precfactor should be at least 2.0.  Changing precfactor to 2.0');
	PDEparam.precfactor = 2.0;
end
if ~rval
    return;
end

myfontsize=16;
mysmallfontsize=14;
%points = 144*3; %spacing between labels on contours - made so that only one label appears per line

% This code follows the idea of Chernoff & Petkau, and of Brezzi and Lai,
% papers, except that we use a trinomial tree rather than binomial tree.

% Start setting up computational variables
if PDEparam.finiteT
    sinit = sEND;      % start the calculation of the free boundary at values beyond sEND (smaller than sEND) in reverse time scale
else
    sinit = 0.98*sEND;      % start the calculation of the free boundary at values beyond sEND (smaller than sEND) in reverse time scale
end
if isa(approxmeth, 'function_handle')
%    if ~isDisc % for discounted case, want to have precfactor>1 grid points above 0 where expected reward for continuing is at least 0
%        dw = sinit^2/(20/9)/PDEparam.precfactor;  % need to divide be something bigger than 2 in order to get a valid recursion when C&F model is used
%    else       % for discounted case, want at least precfactor>1 grid points below the asymptotic approx for boundary for small sinit
%        dw = (min(0.2,sinit)/sqrt(2))/PDEparam.precfactor/(10/9);      % need to divide be something bigger than 2 in order to get a valid recursion when C&F model is used
%    end
    dw = approxmeth(sinit,PDEscale,PDEparam) / PDEparam.precfactor / (1.5);    % the extra 10/9 is a fudge factor for hoping for numerical error improvement
else
    dw = approxmeth(1);
end
ds = dw^2 * 2 / 3;          % in trinomial tree, equal probs of going up, straight or down implies this equation to get correct variance of reverse time brownian motion
%SAVEEVERY=280;        % every how many iterations do we keep the PDE values, for storage in matrix and files?
%NUMSAVESPERITER=100;  % number of times we save the PDE values, for storage in matrix and files?
SAVEEVERY=350;        % every how many iterations do we keep the PDE values, for storage in matrix and files?
NUMSAVESPERITER=80;  % number of times we save the PDE values, for storage in matrix and files?
Numds=SAVEEVERY*NUMSAVESPERITER;   % number of time steps per grid, % 8000 or 16000 for example, 
FRACTOKEEP = 0.9; % use this to remember the state at some earlier time: 0.9 or 085 should get of ripples in many cases, make sure it is smaller than
if (FRACTOKEEP > (1-2/SAVEEVERY)) || (FRACTOKEEP < 0.5)
    warning('keep FRACTOKEEP between 0.5 and 1 - 2/SAVEEVERY, please insure SAVEEVERY big enough');
end
myeps = 10^-10; %min(1,beta)*min(10^-12, ds/1000);    % a small value, for use in checking if things are close to 0 or not


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHUNK1: Set up the Grid

% now figure out what is the right number of dw values which are needed
% along the way - knowing that dw and ds will be scaled up, and the
% boundaries are roughly approximated using the asymptotic
% approximations...
THRESHSCALE = 100;          % set to be a bit bigger than the threshold value for 's' in the scaled approx for optimal boundary
numchunks = 1 + ceil( log2( 1.1*max(1,min(THRESHSCALE,s0)-sinit)/(Numds*ds) ) ); % guess as to how many iterations are needed to generate all the files
spowvec = [0 4.^(0:max(0,(numchunks-1)))]; % compute some vectors to say roughly how many times the initial 'ds' are the s-scale steps in the i=1,...,numchunks iterations
wpowvec = [1 2.^(0:max(0,(numchunks-1)))]; % compute some vectors to say roughly how many times the initial 'dw' are the w-scale steps in the i=1,...,numchunks iterations
smaxvec = sinit+Numds*ds*cumsum(spowvec); 
% Right now we are only checking for figure fitting within boundary. Maybe
% we should also be checking for number of grid points within boundary at
% each value of the boundary... might need to shrink dw again...
if isa(approxmeth, 'function_handle')
    wmaxvec=approxmeth(smaxvec,PDEscale,PDEparam);      % approximate the maximum value of w for which the grid must account
    dwvec = dw*wpowvec;
    %dsvec = ds*spowvec;
    ceilfactor=PDEparam.ceilfactor; % multiple by 4 because of granularity and fact that we are not checking all points for (s, b(s)), just those at grid size changes
    ratiovec = wmaxvec ./ dwvec;    % figure roughly how many dw from 0 to upper boundary...
    bigw = ceil(ceilfactor*max(ratiovec)); % hit that figure with a multiplier in order to make sure there is a margin of error, especially below boundary
else  % if we don't have a function to compute, approximately, the upper stopping boundary, then take it from user specified input
    bigw = max(2*PDEparam.precfactor,ceil(approxmeth(2)));
end

% Now, can initialize the vector of w values for the recursion, and compute
% how big it is.
wvec = dw*(-bigw:bigw);
[wvsize] = length(wvec);

hifrac = 0; lowfrac = 1; % These will be used in diagnostics to see if ceilfactor is set up well or not.
% if there is no lower boundary, then it is fine if lowfrac is 1 at the end
% of a run. however, when discount rate is 0, then hifrac and lowfrac
% should both end up near .5 to .9. if they are a lot smaller than that at
% the end of a run, then ceilfactor can be made smaller for the case of no
% discounting - this will help the code run faster.

% Set up the initial conditions for the recursions
scur = sinit;

% vectors for value function
[Cinitvec, sincrem] = rewardfunc(wvec,scur,PDEscale, PDEparam); % Approximation of value to go, given stopping at the discretized time sinit.
Cin=Cinitvec;    % initialize terminal condition assuming immediate stopping.
Cout=Cin;           % only to initialize
% vectors for EN = expected number of samples
ENCin = sincrem; %0*Cinitvec; % a priori, assume 0 more samples in expectation at terminal time
ENCout = ENCin; % a priori, assume 0 more samples in expectation at terminal time
% structures for PCS
EPCSCout = normcdf(abs(wvec)./sqrt(scur),0,1); % probability of correct selection if choice made at time scur
EPCSCin = EPCSCout; 

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
if PDEparam.DoPlot
    middlecin1 = Cinitvec((bigw-4):(bigw+5))-max(0,wvec((bigw-4):(bigw+5)));
    disp(middlecin1);
end
minindx=1; counter = 1;
while minindx==1
    % do the iteration
    sout = scur+ds;
    Cinitvec = termreward(wvec,sout,PDEscale, PDEparam);  % compute reward from stopping
    if ~isDisc % for undiscounted rewards, do the update, noting that the top and bottom can be 'tied' to be 'stopped' 
        Cout(2:(wvsize-1)) = ( Cin(1:(wvsize-2)) + Cin(2:(wvsize-1)) + Cin(3:wvsize)) / 3  - (1/scur-1/sout);
        if ~isOffline   % if there is online learning, add a term to that effect
            Cout(2:(wvsize-1)) = Cout(2:(wvsize-1)) + (1/scur-1/sout) * wvec(2:(wvsize-1));
        end
        Cout(1) = Cinitvec(1);              % following conditions in sure stopping at top and at bottom of vector
        Cout(wvsize) = Cinitvec(wvsize);
    else  % with discounting, it might be possible to have no lower stopping region, therefore we need to handle the update for Cout(1) in special way
        % we use a surrogate for the reward of doing dw below the stopping boundary... 
        Cout(2:(wvsize-1)) =  exp(- (1/scur-1/sout)) * ( Cin(1:(wvsize-2)) + Cin(2:(wvsize-1)) + Cin(3:wvsize)) / 3 ;
        if ~isOffline   % if there is online learning, add a term to that effect
            Cout(2:(wvsize-1)) = Cout(2:(wvsize-1)) + (1-exp(-(1/scur-1/sout))) * wvec(2:(wvsize-1));
        end
        Cout(1) = exp(- (1/scur-1/sout)) * (termreward(min(wvec)-dw,sout,PDEscale, PDEparam) + Cin(1) + Cin(2) ) / 3 ;   % following conditions in sure stopping at top and at bottom of vector
        Cout(wvsize) = Cinitvec(wvsize);
    end
    Cout=max(Cinitvec,Cout);
    Cin=Cout;                           % iterate to obtain new initial condition

    stopped = ( abs(Cout - Cinitvec) < myeps ); % check if we have stopped or not.
    [~, minindx] = min(stopped);             % get index to first grid point above lower stopping boundary
    [~, maxindx] = min(flipdim(stopped,2));  % get index to first grid point below upper stopping boundary
    maxindx = wvsize + 1 - maxindx; %(scur < sEND - 2*ds)   % fix the value of minindx
    hifrac = max(hifrac,maxindx/length(stopped)); % compute some diagnostic stats for how much of the grid is being used by the contin set.
    lowfrac = min(lowfrac,minindx/length(stopped));

    EPCSCout(2:(wvsize-1)) = ( EPCSCin(1:(wvsize-2)) + EPCSCin(2:(wvsize-1)) + EPCSCin(3:wvsize)) / 3; % update PCS
    EPCSCout(stopped) = normcdf(abs(wvec(stopped))./sqrt(sout),0,1);
    EPCSCin = EPCSCout; 
    ENCout(2:(wvsize-1)) = ( ENCin(1:(wvsize-2)) + ENCin(2:(wvsize-1)) + ENCin(3:wvsize)) / 3  + (1/scur-1/sout);
    ENCout(stopped) = 0;
    ENCin = ENCout; 

    if (~PDEparam.finiteT) && (sout < sEND - 2*ds) % unless it is for finite horizon, force a few iterations of this recursion in order to 'prime' the initial conditions and to avoid 'edge' effects from the  initializaiton of the terminal condition
%    if (counter < 20*PDEparam.precfactor)  % force a few iterations of this recursion in order to 'prime' the initial conditions and to avoid 'edge' effects from the initializaiton of the terminal condition
        minindx = 1;
    end % sEND
    
    scur=sout;
    counter = counter + 1;
end 
if PDEparam.DoPlot
    disp(dw);
    disp(ds);
    disp(bigw);
    disp(scur);
    middlecin = Cin((bigw-4):(bigw+5))-max(0,wvec((bigw-4):(bigw+5)))-Cinitvec((bigw-4):(bigw+5))-max(0,wvec((bigw-4):(bigw+5)));
    disp(middlecin);
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
ijk=0;   % index to keep track of the number of blocks (one per rescaling of dw and ds) used in the computations
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
%            Cout(1) = exp(- (1/scur-1/sout)) * (rewardfunc(min(wvec)-dw,sout,PDEscale, PDEparam) + Cin(1) + Cin(2) ) / 3 ;              % following conditions in sure stopping at top and at bottom of vector
            Cout(1) = exp(- (1/scur-1/sout)) * (termreward(min(wvec)-dw,sout,PDEscale, PDEparam) + Cin(1) + Cin(2) ) / 3 ;              % following conditions in sure stopping at top and at bottom of vector
            Cout(wvsize) = Cinitvec(wvsize);
        end
        Cout=max(Cinitvec,Cout);
        Cin=Cout;                           % iterate to obtain new initial condition

        % try to compute the stopping boundaries
        stopped = ( abs(Cout - Cinitvec) < myeps );
        [~, minindx] = min(stopped);             % get index to first grid point above lower stopping boundary
        [~, maxindx] = min(flipdim(stopped,2));  % get index to first grid point below upper stopping boundary
        maxindx = wvsize + 1 - maxindx;         %    [scur minindx maxindx];
        hifrac = max(hifrac,maxindx/length(stopped)); % compute some diagnostic stats for how much of the grid is being used by the contin set.
        lowfrac = min(lowfrac,minindx/length(stopped));

        % update expected number of samples til stopping time: to be
        % debugged
        ENCout(2:(wvsize-1)) = ( ENCin(1:(wvsize-2)) + ENCin(2:(wvsize-1)) + ENCin(3:wvsize)) / 3  + (1/scur-1/sout);
        ENCout(stopped) = 0;
        ENCin = ENCout; 

        % compute PCS: to be debuged. 
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
    mymat = strcat(matdir,fName,int2str(ijk),'.mat');
% should be able to reconstruct Vwsmatrix from Bwsmatrix and PDEscale, PDEparam.
    save(mymat,'Bwsmatrix','ENwsmatrix','EPCSwsmatrix','svec','wvec','upvec','downvec','up1','down1');

    % if diagnostic plots are desired, plot them
    if PDEparam.DoPlot
        disp(dw);
        disp(ds);
        disp(hifrac);
        disp(lowfrac);
        %        figure(180+ijk)
        if ~exist('fignum','var'), fignum = 20; end;
        fignum=fignum+1;figure(fignum);
        hold off
        plot(wvec(:),Cinitvec(:),'-.',wvec(:),Cin(:),'--');set(gca,'FontSize',mysmallfontsize);
        hold on; 
        plot(downvec(ind),0,'x',upvec(ind),0,'x');set(gca,'FontSize',mysmallfontsize);
        % tmp=axis;tmp(4)=1.2*max(max(10*dw,max(up1))/beta,max(upvec)/beta);tmp(3)=0;axis(tmp);
        tmp=axis;tmp(1)=(min(down1)-dw);tmp(2)=(max(up1)+dw);tmp(3)=0;axis(tmp);
        xlabel('wvec','FontSize',myfontsize,'FontName','Times'); ylabel('V(s,w)','FontSize',myfontsize,'FontName','Times')
        title('Value function in (w,s) coords','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigVWS',int2str(ijk),'.eps');
        if figsave 
            print('-deps',mytitle); 
        end	

        if ~exist('fignum','var'), fignum = 20; end;
        fignum=fignum+1;figure(fignum);
        hold off
        semilogy(wvec(:),Cinitvec(:),'-.',wvec(:),Cin(:),'--');set(gca,'FontSize',mysmallfontsize);
        %plot(wvec(:),Cinitvec(:),'-.',wvec(:),Cin(:),'--');set(gca,'FontSize',mysmallfontsize);
        hold on; 
        legend('Cinitvec','Cin');
        semilogy(downvec(ind),10^-3,'x',upvec(ind),10^-3,'x');set(gca,'FontSize',mysmallfontsize);
        %plot(downvec(ind),10^-3,'x',upvec(ind),10^-3,'x');set(gca,'FontSize',mysmallfontsize);
        tmp=axis;tmp(1)=(min(down1)-dw);tmp(2)=(max(up1)+dw);tmp(3)=0;axis(tmp);
        xlabel('wvec','FontSize',myfontsize,'FontName','Times'); ylabel('V(s,w)','FontSize',myfontsize,'FontName','Times')
        title('Value function in (w,s) coords','FontSize',myfontsize,'FontName','Times')
        mytitle = strcat(figdir,fName,'FigLogVWS',int2str(ijk),'.eps');
        if figsave 
            print('-deps',mytitle); 
        end	
        WidthContin(ijk) = maxindx - minindx;
    end  % diagnostic plot routines

    % double check boundaries to see if bias corrected bound exceeds max w
    % value in calculation or not. if there is an error here, then need to
    % set bigw to a larger value
    if max(max(up1),max(up1)) >= max(wvec) % also put -down1?
        warning('error: need bigger initial number of grid points for w direction');
%        return;
    end

    %   now reset the ds and dw values to scale up an order of magnitude
    %FIX: if code to be extended for further orders of magnitude difference
    %between s0 and sEND, then this might be made adaptive: doing the
    %doubling of dw only if there is still enough number of dw up to
    %boundary, might need to increase bigw if this adaptive grid
    %flexibility is added.
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
mymat = strcat(matdir,fName,int2str(ijk),'.mat');
save(mymat,'fName', 'TimeStamp','StartFileVal','EndFileVal','PDEscale','PDEparam','lasts','firsts','lasthelds','myeps','hifrac','lowfrac','WidthContin');

end