function [ evivec, dsvec ] = PDECGApproxValue( wvec, s0, scale )
%PDECGApproxValue: Find the best on can do in terms of maximizing expected reward
% for the undiscounted standardises stopping problem in the standardized
% (w,s) reverse time scale, for the stopping problem:
%   max E[ exp(- (1/S - 1/s0)) W_S | w_{s_0} = wvec, s_0 ]
% where the maximization is with respect to stopping times S, on compact
% interval starting at s0 and running backward to 0, and where S can not
% depend on the value of W (it must be selected in advance). Returns the
% maximizing reward in evivec for each w in wvec, and the change in S down
% from s0 in dsvec, assuming that scale has the scale for the original
% problem (with coefficient gamma as a field) or with gamma passed as
% scale.
%
% Note that this is NOT the KG* of Frazier and Powell adapted to discounted
% rewards. This function seeks to maximise expected reward, NOT find the 
% sampling plan which maximizes reward per sample.
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
    NUMCHECKS = 200; % FIX: Can probabably find a way to speed this check, but this is reasonable proxy for the moment.
        % try to find KG* type 'best' lookahead value in (w,s) scale when
        % w0=0 and s=s0
    svec = s0*(1-(1:NUMCHECKS)/(3*NUMCHECKS));
    sincrem = s0-svec(svec>0 & svec < s0);             % time elapse from scur to valid values in svec    
    for i=1:length(sincrem)           % check the lookaheads over that interval
        %erew = (max(0,wvec+dincrem(i)) + max(0,wvec) + max(0,wvec-dincrem(i)))/3; % reward making decision in a bit of time
        %evitmp = exp(1/s0 - 1/svec(i)) * erew; % this is in spirit of
        evitmp = exp(1/s0 - 1/(s0-sincrem(i))) * sqrt(sincrem(i)) * PsiNorm(-wvec/sqrt(sincrem(i)));
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
    HIVAL = 16;      % check lookahead up to 2^HIVAL steps aheads
    NUMCHECKS = 200; % check 20 values of increment
    tstincrem = 2.^(LOWVAL:((HIVAL-LOWVAL)/NUMCHECKS):HIVAL);
    svec = s0 * 1 ./ ( 1 + tstincrem*gamma*s0);
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

