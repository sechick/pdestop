function [ evivec ] = PDECFKGs( wvec, s0, svec )
%PDECFKGs:  Solves something similar to the KG* in the (w, s) standardized
%reverse time scale, for the stopping problem:
%   max E[ - (1/S - 1/s0) + max(0, W_S) | w_{s_0} = wvec, s_0 ]
% for each value of the standardized mean in wvec, and where stopping times
% are tested for S in svec \cup s0(for values of svec in (0, s0]
% If svec is passed as empty string, then a default vector of values of
% svec is created.
%
% This problem arises on page 568 of C&F (2012). The function values are
% useful for approximating the reward function of the free boundary
% problem at the time sEND (smallest value of -s solved for in the
% diffusion) and is a 'better' approximation than the old default value of
% max(wvec,0). The philosophy is one of approximating the reward function
% with a KG*-like test of multiple one-step lookahead durations whose
% values are mapped to svec.

% At worst case, assume that one can stop immediately and get the best of 0
% and the value of wvec. This is the usual terminal reward for stopping.
evivec = max(wvec,0);   % unless otherwise specified, best you

% do some error checking
if s0 <= 0
    warning('s0 should be strictly greater than 0');
else
    NUMSVECDEFAULT = 100;
    if nargin < 3 % if svec is not passed, then implement a 'default' set of values for it
        % try to find KG* type 'best' lookahead value in (w,s) scale when
        % w0=0 and s=s0
        dwtest = min(1,s0)^2/2/10;
        dstest = dwtest^2*2/3;
        svec = s0 - dstest*(1:NUMSVECDEFAULT)/10;
    end
    svec2 = svec(svec>0 & svec < s0);       % pull out the values of svec in open interval (0, s0)
    for i=1:length(svec2)
        evitmp = -(1/svec2(i) - 1/s0) + sqrt(s0 - svec2(i)) * PsiNorm(-wvec/sqrt(s0 - svec2(i)));
        evivec = max(evivec,evitmp);
    end
end
    
end

