function [rval] = PsiNormRepsRoot(y,sig,t,reps,value)
% (c) 2004 Stephen E. Chick, all rights reserved
% This file is for input to matlab, and does calculations
% to support the 'selecting the best system' paper in the
% bayesian environment.
%
% Code is 'as is' and no guarantees for correctness.
%
% INPUTS:
% y = tested value of y_t
% sig = standard deviation of samples
% t = effective number of samples
% value = target EOC bound
%
% OUTPUTS:
% the difference in actual and targeted value for EOC
%
% Example use: 
%   fzero(@PsiNormRoot,1,optimset('TolX',1e-8),1,2,.01)
%   returns the value of y such that the target EOC is used.
%   that means that y/t (the value returned for y / t) is the stopping
%   boundary for the EOC-based stopping rule
%
%[y sig t value]
%mystder = sig * reps / sqrt( t * (t+reps));
mystder = sig * sqrt(reps) / sqrt( t * (t+reps));
mymean = y/t;
myEVI = mystder * PsiNorm( mymean / mystder );
rval = myEVI - value;