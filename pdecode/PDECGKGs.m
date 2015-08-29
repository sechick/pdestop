function [ evivec, dsvec ] = PDECGKGs( wvec, s0, scale )
%PDECFKGs:  Solves something similar to the KG* in the (w, s) standardized
%reverse time scale, for the stopping problem:
%   max E[ exp(- (1/S - 1/s0)) W_S | w_{s_0} = wvec, s_0 ]
% for each value of the standardized mean in wvec, and where stopping times
% are tested for S in svec \cup s0 (for values of svec in (0, s0]
% If scale is not passed as an argument, then a default vector of values of
% svec is created. If scale is passed as a real, then it is interpreted as
% gamma, the rescacing parameter s = 1/(gamma n_t), for the -s scale
% diffusion. If scale is passed as a structure, then the .gamma field is
% used for the value of gamma.
%
% This problem arises on page ec7 in equation (ec.13) of C&G (2009). 
% The function values are
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
        %erew = (max(0,wvec+dincrem(i)) + max(0,wvec) + max(0,wvec-dincrem(i)))/3; % reward making decision in a bit of time
        %evitmp = exp(1/s0 - 1/svec(i)) * erew; % this is in spirit of
        evitmp = exp(1/s0 - 1/svec(i)) * sqrt(sincrem(i)) * PsiNorm(-wvec/sqrt(sincrem(i)));
        dsvec(evitmp > evivec) = sincrem(i);
        evivec = max(evivec,evitmp);
    end
else  % third argument was passed - try to pick out scaling parameter gamma from it then compute KG*
    if isfield(scale,'gamma') % need gamma in order to figure out what ds for one sample
        gamma = scale.gamma;
    else
        gamma = scale;
    end
    
    LOWVAL = 0;     % require at least 2^0 = 1 step look ahead
    HIVAL = 7;      % check lookahead up to 2^HIVAL steps aheads
    NUMCHECKS = 20; % check 20 values of increment
    tstincrem = 2.^(LOWVAL:((HIVAL-LOWVAL)/NUMCHECKS):HIVAL);
    svec = s0 * 1 ./ ( 1 + tstincrem*gamma*s0)
    sincrem = s0-svec(svec>0 & svec < s0) ;            % time elapse from scur to valid values in svec    

    for i=1:length(sincrem)           % check the lookaheads over that interval
%        erew = (max(0,wvec+dincrem(i)) + max(0,wvec) + max(0,wvec-dincrem(i)))/3; % reward making decision in a bit of time
%        evitmp = exp(1/s0 - 1/svec(i)) * erew;
        evitmp = exp(1/s0 - 1/(s0-sincrem(i))) * sqrt(sincrem(i)) * PsiNorm(-wvec/sqrt(sincrem(i)));
        dsvec(evitmp > evivec) = sincrem(i) ;
        evivec = max(evivec,evitmp);
    end
    % don't yet have bounds for KG* approx for discounted case which would
    % be analogous to the bounds for undiscounted in frazier and powell
    % 2010
end
    
end

