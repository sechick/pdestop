function [ Vvec, Bwsval, ENvec, PCSvec ] = PDEGetVals(cfSoln,wvec,sval) 
%PDEGetVals: Extrapolates sPDE solution based on computations loaded in
%cfSoln, at vector of w values in wvec, at time s=sval. That is, this is
%the expected reward computed by the PDE solver, when the unkown mean upon
%starting, in the standardized (w, s) scale, is Normal(wvec, sval)
%
% requires sval to be SCALAR!
%
% OUTPUTS: Vvec contains the value function, Bwsval contains the expected
% value of sequential sampling (should be 0 outside of the continuation
% region), ENvec contains a scaled version of the expected number of
% samples (divide by beta to get right number of sample), PCSvec contains
% the probability of correct selection, in expectation, given wvec, sval.
%
% For purposes here, continuation region is approximated by the upper and
% lower boundary computed by the finite difference method, not the
% asymptotic approximations
%
% it is presumed that cfSoln was loaded by PDESolnLoad(), and that function
% in turn presumes that PDESolnCompute had been used to compute the PDE
% solution
%
% (c) 2015, S Chick, all rights reserved. Code is provided 'as is' with no
% warranty of correctness.

if (length(sval) > 1)
    warning('unpredictable behavior in PDEGetVals if sval is not a scalar');
end
% Do a bit of data validation
% load in data from solution which will be required to compute the return
% values
stopfunc = cfSoln.Header.PDEparam.termrewardfunc;
approxfunc = cfSoln.Header.PDEparam.approxvaluefunc;
firsts = cfSoln.Header.firsts;
lasts = cfSoln.Header.lasts;
biggests = lasts(length(lasts));

% OUTPUT VALUE INITIALIZATION
% initialize Vvec, the output with value of continuing, to be the value of
% stopping immediately. 
Vvec = stopfunc(wvec(:),sval,cfSoln.Header.PDEscale,cfSoln.Header.PDEparam);
% initialize value of information of continuing to sample to 0, and PCS to
% be the PCS given stopping immediately... 
Bwsval = 0*Vvec;
PCSvec = normpdf(abs(wvec(:))/sqrt(sval),0,1);
% initialize the other outputs to the empty vector: their values will be
% filled in if the can be filled in below
ENvec = [];
% The values of these outputs will be updated below if it is possible to do so

% OUTPUT VALUE UPDATE IF POSSIBLE WITH SAVED FILES
% Find which block to use to estimate function value at time sval
if sval < firsts(1)  % sval was below the range for which data was originally computed
    % therefore, use approximation to value function provided during
    % computation time 
    Vvec = approxfunc(wvec(:),sval,cfSoln.Header.PDEscale,cfSoln.Header.PDEparam);
    Bwsval = Vvec-stopfunc(wvec(:),sval,cfSoln.Header.PDEscale,cfSoln.Header.PDEparam);
    % FIX: Might be able to do a smarter approximation by looking at
    % solution at smallest value of s which had computations.
    warning( sprintf('sval = %f is below range of s, (%f, %f), computed: using approx func',sval,firsts(1),biggests));
elseif sval > biggests % sval is above the range for which data was originally computed
    Vvec = approxfunc(wvec(:),sval,cfSoln.Header.PDEscale,cfSoln.Header.PDEparam);
    Bwsval = Vvec-stopfunc(wvec(:),sval,cfSoln.Header.PDEscale,cfSoln.Header.PDEparam);
    warning( sprintf('sval = %f is above range of s, (%f, %f), computed: using approx func',sval,firsts(1),biggests));
else
    isabove = (sval > lasts);
    [~, ind] = min(isabove);      % here, we know minimum will be 0 becase sval <= biggests, so b is index of file to use
    wv = cfSoln.Data(ind).wvec;
    svec = cfSoln.Data(ind).svec;
    Bwsmatrix = cfSoln.Data(ind).Bwsmatrix;
    upper = interp1( svec, cfSoln.Data(ind).upvec, sval );
    lower = interp1( svec, cfSoln.Data(ind).downvec, sval );
    Bwsval=interp2(svec,wv,Bwsmatrix,sval,wvec(:));
    if min(Bwsval) < 0
        warning('interpolated value of continuing to sample computed to be negative')
    end
    Bwsval((wvec < upper) & (wvec > lower)) = abs(Bwsval((wvec < upper) & (wvec > lower)));
    Vvec((wvec < upper) & (wvec > lower)) = Vvec((wvec < upper) & (wvec > lower)) + Bwsval((wvec < upper) & (wvec > lower));      % updated value for Vvec is value of continuing to sample plus value of stpping now
    ENvec = interp2(svec,wv,cfSoln.Data(ind).ENwsmatrix,sval,wvec(:));
    ENvec(((wvec >= upper) & (wvec <= lower))) = 0.0;
    PCSvec = interp2(svec,wv,cfSoln.Data(ind).EPCSwsmatrix,sval,wvec(:));
end


end

