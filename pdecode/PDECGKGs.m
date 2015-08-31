function [ kgfactor, dsvec ] = PDECGKGs( wvec, s0, scale )
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
% Kgfactor is returned somewhat analogous to what has been done for offline
% discounted kg factor by ryzhov, powell, frazier, for infinite horizon. It
% assumes one retains mu for some period of time, then one observes the
% output, and then one picks the best, suitably discounted. This is done in
% the (w,s) reverse time coordinate system (which was not treated by RPF).
%
% The value of the optimal one-step lookahead, in -s coordinates, is returned in
% dsvec. This would convert to a lookahead of 1/(gamma ds) samples in t
% coordinates.
%
% At worst case, assume that one can stop immediately and get the best of 0
% and the value of wvec. This is the usual terminal reward for stopping.

FIX

kgfactor = 0*wvec; % initialize the kgfactor to be the value of sticking with the best current alternative (vector of values checked is in wvec)
%kgfactor = max(0,wvec); % initialize the kgfactor to be the value of sticking with the best current alternative (vector of values checked is in wvec)
dsvec = 0*kgfactor;    % initialize the increment for the optimal sampling plan to be 0: it will increase above 0 in the continuation region in calculations below.

if s0 <= 0 % need a valid value of s in order to do computations.
    warning('s0 should be strictly greater than 0');
    return;
elseif nargin < 3   % no parameter passed in for scale: try to flood
    NUMCHECKS = 300; % FIX: Can probabably find a way to speed this check, but this is reasonable proxy for the moment.
        % try to find KG* type 'best' lookahead value in (w,s) scale when
        % w0=0 and s=s0
    svec = s0*(1-(1:NUMCHECKS)/(3*NUMCHECKS));
    sincrem = s0-svec(svec>0 & svec < s0);             % time elapse from scur to valid values in svec    
    for i=1:length(sincrem)           % check the lookaheads over that interval
%        kgtmp = (1 - exp(-sincrem(i)))* wvec + exp(-sincrem(i)) * sqrt(sincrem(i)) * PsiNorm(-wvec/sqrt(sincrem(i)));
        kgtmp = wvec - exp(-sincrem(i)) * sqrt(sincrem(i)) * PsiNorm(-wvec/sqrt(sincrem(i)));
        dsvec(kgtmp > kgfactor) = sincrem(i);
        kgfactor = max(kgfactor,kgtmp);
    end
else  % third argument was passed - try to pick out scaling parameter gamma from it then compute KG*
    if isfield(scale,'gamma') % need gamma in order to figure out what ds for one sample
        gamma = scale.gamma;
    else
        gamma = scale;
    end
    
    LOWVAL = 0;     % require at least 2^0 = 1 step look ahead
    HIVAL = 12;      % check lookahead up to 2^HIVAL steps aheads
    NUMCHECKS = 300; % check 20 values of increment
    tstincrem = 2.^(LOWVAL:((HIVAL-LOWVAL)/NUMCHECKS):HIVAL);
    svec = s0 * 1 ./ ( 1 + tstincrem*gamma*s0)
    sincrem = s0-svec(svec>0 & svec < s0) ;            % time elapse from scur to valid values in svec    

    for i=1:length(sincrem)           % check the lookaheads over that interval
%        kgtmp = (1 - exp(-sincrem(i)))* wvec + exp(-sincrem(i)) * sqrt(sincrem(i)) * PsiNorm(-wvec/sqrt(sincrem(i)));
        kgtmp = wvec - exp(-sincrem(i)) * sqrt(sincrem(i)) * PsiNorm(-wvec/sqrt(sincrem(i)));
        dsvec(kgtmp > kgfactor) = sincrem(i) ;
        kgfactor = max(kgfactor,kgtmp);
    end
    % don't yet have bounds for KG* approx for discounted case which would
    % be analogous to the bounds for undiscounted in frazier and powell
    % 2010
end
    
end

