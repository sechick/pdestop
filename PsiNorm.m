function [rval] = PsiNorm(zval)
% (c) 2004 Stephen E. Chick, all rights reserved
% This file is for input to matlab, and does calculations
% to support the 'selecting the best system' paper in the
% bayesian environment.
%
% Code is 'as is' and no guarantees for correctness.
%
% INPUTS:
% zval : test statistic values (vector)
%
% OUTPUTS:
% rval : vector of outputs
rval = normpdf(zval) - zval .* normcdf( - zval);
rval = max(rval,zeros(size(rval)));
