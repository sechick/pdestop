function [ evivec, dsvec ] = PDECFApproxValue( wvec, s0, scale )
%PDECFKGs: Find the best on can do in terms of maximizing expected reward
% for the undiscounted standardises stopping problem in the standardized
% (w,s) reverse time scale, for the stopping problem:
%   max E[ - (1/S - 1/s0) + max(0, W_S) | w_{s_0} = wvec, s_0 ]
% where the maximization is with respect to stopping times S, on compact
% interval starting at s0 and running backward to 0, and where S can not
% depend on the value of W (it must be selected in advance). Returns the
% maximizing reward in evivec for each w in wvec, and the change in S down
% from s0 in dsvec, assuming that scale has the scale for the original
% problem (with coefficient gamma as a field) or with gamma passed as
% scale.
%
% Note that this is NOT the KG* of Frazier and Powell: it seeks to maximise
% expected reward, NOT find the sampling plan which maximizes reward per
% sample.
%
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
    svec = s0 * 1 ./ ( 1 + tstincrem*gamma*s0);
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
end


end
