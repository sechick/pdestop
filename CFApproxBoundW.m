function [wval] = CFApproxBoundW(sval)
% CFApproxBoundW: Computes approximate upper boundary for optimal stopping
% time for Bayesian bandit problem, as in Chick & Frazier Mgmt Sci 2012,
% for the case of positive sampling costs and zero discount rate.
%
% Computations done in standardize reverse time brownian motion
% coordinates. That is, in the (w,s) scale, where the point (w,s) refers to
% a prior distribution for an unknown mean which is normal(w, s).
%
% INPUTS:
%   sval: vector of values for evaluating the approximate boundary
% OUTPUTS:
%   wval: vector of bountary values corresponding to the inputs
%
% This function is based on equation (16) of Chick & Frazier (Mgmt Sci
% 2012), with an update for values of sval > 40 to take advantage of
% asymptotic approximations of Chernoff (Ann Math Statist, 1965, paper on
% sequential sampling III).
%
% Code is 'as is' and no guarantees for correctness.
%
% (c) 2009 Stephen E. Chick, all rights reserved
%   Updated 2015 Aug 14, to bring in Chernoff Sequential test III result
%   for boundary for very large sval.

    wval = sval;
    wval(sval <= 1) = 0.233*sval(sval<=1).^2;
    wval(sval>1 & sval<=3) = 0.00537*sval(sval>1 & sval<=3).^4 - 0.06906*sval(sval>1 & sval<=3).^3 + ...
        0.3167*sval(sval>1 & sval<=3).^2 - 0.02326*sval(sval>1 & sval<=3);
    wval(sval>3 & sval<=40) = 0.705*sqrt(sval(sval>3 & sval<=40)).*log(sval(sval>3 & sval<=40));
    wval(sval > 40) = sqrt( sval(sval > 40) .* (3*log(sval(sval > 40))) - log(8*3.14159) - 2./log(sval(sval > 40)) - 170);

end     