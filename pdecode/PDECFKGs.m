function [ evivec ] = PDECFKGs( wvec, s0, scale )
%PDECFKGs:  Solves something similar to the KG* in the (w, s) standardized
%reverse time scale, for the stopping problem:
%   max E[ - (1/S - 1/s0) + max(0, W_S) | w_{s_0} = wvec, s_0 ]
% for each value of the standardized mean in wvec, and where stopping times
% are tested for S in svec \cup s0(for values of svec in (0, s0]
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

% At worst case, assume that one can stop immediately and get the best of 0
% and the value of wvec. This is the usual terminal reward for stopping.
evivec = max(wvec,0);   % lower bound on best you can do is to stop immediately

if s0 <= 0 % need a valid value of s in order to do computations.
    warning('s0 should be strictly greater than 0');
    return;
elseif nargin < 3   % no parameter passed in for scale: try to flood
    NUMCHECKS = 300; % FIX: Can probabably find a way to speed this check, but this is reasonable proxy for the moment.
        % try to find KG* type 'best' lookahead value in (w,s) scale when
        % w0=0 and s=s0
    svec = s0*(1-(1:NUMCHECKS)/(4.5*NUMCHECKS));
    sincrem = s0-svec(svec>0 & svec < s0);             % time elapse from scur to valid values in svec    
    dincrem = sqrt(3*sincrem/2);    % implied reward from diffusion over that increment
    for i=1:length(dincrem)           % check the lookaheads over that interval
        evitmp = -(1/svec(i) - 1/s0) + sqrt(sincrem(i)) * PsiNorm(-wvec/sqrt(sincrem(i)));
        evivec = max(evivec,evitmp);
    end
else    
    if isfield(scale,'gamma') % need gamma in order to figure out what ds for one sample
        gamma = scale.gamma;
    else
        gamma = scale;
    end
    svec = s0*(1-gamma*2.^(0:.1:6)); % check kg_x where x = 2.^(0:.2:8)
    sincrem = s0-svec(svec>0 & svec < s0) ;            % time elapse from scur to valid values in svec    
    dincrem = sqrt(3*sincrem/2);    % implied reward from diffusion over that increment
    for i=1:length(dincrem)           % check the lookaheads over that interval
        evitmp = -(1./(svec(i)) - 1/s0) + sqrt(sincrem(i)) .* PsiNorm(-wvec/sqrt(sincrem(i)));
        evivec = max(evivec,evitmp);
    end
    
    % now check bounds from Frazier and Powell 2010
    kgstarterm = wvec.^2 / s0; 
    sfactor = (kgstarterm - 1 + sqrt(kgstarterm.^2 + 6*kgstarterm + 1)) / 4; % m-underbar from prop 4, frazier & powell, 2010
    sincrembase = max( gamma / (1 + gamma), sfactor ./ ( (1 + sfactor) * s0)); % this is the increment for kgstar for C&F. 
	evitmp = -(1./(s0-sincrembase) - 1/s0) + sqrt(sincrembase) .* PsiNorm(-wvec./sqrt(sincrembase));
    evivec = max(evivec,evitmp);

    sfactor = (kgstarterm - 1 + sqrt(kgstarterm.^2 + 14*kgstarterm + 1)) / 4; % m-underbar from prop 4, frazier & powell, 2010
    sincrembase = max( gamma / (1 + gamma), sfactor ./ ( (1 + sfactor) * s0)); % this is the increment for kgstar for C&F. 
	evitmp = -(1./(s0-sincrembase) - 1/s0) + sqrt(sincrembase) .* PsiNorm(-wvec./sqrt(sincrembase));
    evivec = max(evivec,evitmp);
    
    sfactor = (kgstarterm - 1 + sqrt(kgstarterm.^2 + 10*kgstarterm + 1)) / 4; % try middle of bounds, frazier & powell, 2010
    sincrembase = max( gamma / (1 + gamma), sfactor ./ ( (1 + sfactor) * s0)); % this is the increment for kgstar for C&F. 
	evitmp = -(1./(s0-sincrembase) - 1/s0) + sqrt(sincrembase) .* PsiNorm(-wvec./sqrt(sincrembase));
    evivec = max(evivec,evitmp);
end

end
