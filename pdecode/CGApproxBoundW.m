function [wval] = CGApproxBoundW(sval)
% CFApproxBoundW: Computes approximate upper boundary for optimal stopping
% time for Bayesian bandit problem, as in Chick & Gans Mgmt Sci 2009,
% for the case of zero sampling costs and positive discount rate.
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
% This function is based on equation (EC.23) of Chick & Gans (Mgmt Sci
% 2009), with an update for intermediate values of sval to improve
% numerical accuracy.
%
% Code is 'as is' and no guarantees for correctness.
%
% (c) 2009 Stephen E. Chick, all rights reserved
%   Updated 2015 May 10, to improve approximation for intermediate values
%   of sval, assessed with numerical tests during work on Chick, Forster,
%   Pertile (delayed sequential sampling project).

    wval = sval;
    wval(sval <= (1/7)) = sval(sval<= (1/7)) / sqrt(2);
% replacing formula for middle patch with a slightly different version
% which matches the endpoints better (but might call for slightly earlier
% stoppping). If the formula is changed further, please also make the
% change in CGApproxBoundWinv.m
%    wval(sval>(1/7) & sval<=100) = exp( -0.02645*log( sval(sval>(1/7) & sval<=100) ).^2 + ...
%        0.89106* log( sval(sval>(1/7) & sval<=100) ) - 0.4873 );
    wval(sval>(1/7) & sval<=100) = exp( -0.0275*log( sval(sval>(1/7) & sval<=100) ).^2 + ...
        0.8797 * log( sval(sval>(1/7) & sval<=100) ) - 0.5024 );
    wval(sval > 100) = sqrt(sval(sval > 100)) .* sqrt( 2* log(sval(sval > 100)) - ...
        log(log(sval(sval > 100))) - log(16*pi));

     
end    