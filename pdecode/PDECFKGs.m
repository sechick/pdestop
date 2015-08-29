function [ evivec, dsvec ] = PDECFKGs( wvec, s0, scale )
%PDECFKGs:  Solves something similar to the KG* in the (w, s) standardized
%reverse time scale, for the stopping problem:
%   max E[ - (1/S - 1/s0) + max(0, W_S) | w_{s_0} = wvec, s_0 ]
% for each value of the standardized mean in wvec, and where stopping times
% are tested for S in svec \cup s0 (for values of svec in (0, s0]
% If scale is not passed as an argument, then a default vector of values of
% svec is created. If scale is passed as a real, then it is interpreted as
% gamma, the rescacing parameter s = 1/(gamma n_t), for the -s scale
% diffusion. If scale is passed as a structure, then the .gamma field is
% used for the value of gamma.
%
% This problem arises on page 568 of C&F (2012). The function values are
% useful for approximating the reward function of the free boundary
% problem at the time sEND (smallest value of -s solved for in the
% diffusion) and is a 'better' approximation than the old default value of
% max(wvec,0). The philosophy is one of approximating the reward function
% with a KG*-like test of multiple one-step lookahead.
%
% The value of the optimal one-step lookahead, in -s coordinates, is returned in
% dsvec. This would convert to a lookahead of 1/(gamma ds) samples in t
% coordinates
%
% At worst case, assume that one can stop immediately and get the best of 0
% and the value of wvec. This is the usual terminal reward for stopping.
evivec = max(wvec,0);   % lower bound on best you can do is to stop immediately
dsvec = 0*evivec;

if s0 <= 0 % need a valid value of s in order to do computations.
    warning('s0 should be strictly greater than 0');
    return;
elseif nargin < 3   % no parameter passed in for scale: try to flood
    NUMCHECKS = 100; % FIX: Can probabably find a way to speed this check, but this is reasonable proxy for the moment.
        % try to find KG* type 'best' lookahead value in (w,s) scale when
        % w0=0 and s=s0
    svec = s0*(1-(1:NUMCHECKS)/(4*NUMCHECKS));
    sincrem = s0-svec(svec>0 & svec < s0);             % time elapse from scur to valid values in svec    
    for i=1:length(sincrem)           % check the lookaheads over that interval
        evitmp = -(1/svec(i) - 1/s0) + sqrt(sincrem(i)) * PsiNorm(-wvec/sqrt(sincrem(i)));
        dsvec(evitmp > evivec) = sincrem(i);
        evivec = max(evivec,evitmp);
    end
else    
    if isfield(scale,'gamma') % need gamma in order to figure out what ds for one sample
        gamma = scale.gamma;
    else
        gamma = scale;
    end
%    costvec = 0*wvec;

    LOWVAL = 0;     % require at least 2^0 = 1 step look ahead
    HIVAL = 7;      % check lookahead up to 2^HIVAL steps aheads
    NUMCHECKS = 20; % check 20 values of increment
    tstincrem = 2.^(LOWVAL:((HIVAL-LOWVAL)/NUMCHECKS):HIVAL);
    svec = s0 * 1 ./ ( 1 + tstincrem*gamma*s0)
    sincrem = s0-svec(svec>0 & svec < s0) ;            % time elapse from scur to valid values in svec    
%    figure(1001); hold off; figure(1002); hold off; figure(1003); hold off;
    for i=1:length(sincrem)           % check the lookaheads over that interval
        evitmp = -(1./(s0-sincrem(i)) - 1/s0) + sqrt(sincrem(i)) .* PsiNorm(-wvec/sqrt(sincrem(i)));
        dsvec(evitmp > evivec) = sincrem(i) ;
%        costtmp = -(1/(s0-sincrem(i)) - 1/s0);
%        costvec(evitmp > evivec) = costtmp;
        evivec = max(evivec,evitmp);
%        figure(1001); plot(wvec,sincrem(i)); hold on; plot(wvec,s0);
%        figure(1002); plot(wvec,evivec,'k',wvec,evitmp,'b'); hold on;
%        figure(1003); plot(wvec,costvec); hold on;
    end
    
    % now check bounds from Frazier and Powell 2010
    evivec2 = max(wvec,0);   % lower bound on best you can do is to stop immediately
    dsvec2 = 0*evivec2;   % lower bound on best you can do is to stop immediately
    NUMKGCHECKS = 12;
    CHECKVEC=6:(4/NUMKGCHECKS):10;
    kgstarterm = wvec.^2 / s0; 
%    figure(10001); hold off; figure(10002); hold off;
    for i=1:length(CHECKVEC)
        sfactor = ( 1 /(4*gamma*s0)  ) * (kgstarterm - 1 + sqrt(kgstarterm.^2 + CHECKVEC(i)*kgstarterm + 1)); % m-underbar from prop 4, frazier & powell, 2010
        sincrembase = s0 * max( gamma*s0 / (1 + gamma*s0) , sfactor*gamma*s0 ./ ( 1 + sfactor*gamma*s0) ); % this is the increment for kgstar for C&F. 
        evitmp = -(1./(s0-sincrembase) - 1/s0) + sqrt(sincrembase) .* PsiNorm(-wvec./sqrt(sincrembase));
        dsvec2(evitmp > evivec2) = sincrembase(evitmp > evivec2);%(evitmp > evivec2);
        evivec2 = max(evivec2,evitmp);
%        figure(10001); plot(wvec,sincrembase); hold on; plot(wvec,s0);
%        figure(10002); plot(wvec,evivec2,'k',wvec,evitmp,'b'); hold on;
%        figure(10003); plot(wvec,costvec); hold on;    end
    evivec-evivec2;
%    figure(2000)
%    plot(wvec,evivec,'-',wvec,evivec2,'--')
%    figure(2001)
%    plot(wvec,dsvec,'-',wvec,dsvec2,'--')

    evivec = max(evivec,evivec2);
    dsvec = max(dsvec,dsvec2);

end


end
