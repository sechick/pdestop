function [ Vvec, Bwsval, ENvec, PCSvec, errstat ] = PDEGetVals(cfSoln,wvec,sval,extraflag) 
%PDEGetVals: 
% Description and INPUTS:
% Extrapolates sPDE solution based on computations loaded in cfSoln (from
% PDE solver) at vector of w values in wvec, at time s=sval. That is, the
% the expected reward computed by the PDE solver, when the unkown mean upon
% starting, in the standardized (w, s) scale, is Normal(wvec, sval).
%
% requires sval to be SCALAR! so only wvec values at a single value of sval
% can be returned per call.
%
% extraflag: optional argument, defaults to false. If set to true, then
% ENvec and PCSvec are also returned. Otherwise they are returned as the
% empty vector.
%
% OUTPUTS: Vvec contains the value function, Bwsval contains the expected
% value of sequential sampling (should be 0 outside of the continuation
% region). If extraflag is true, then ENvec contains a scaled version of 
% the expected number of samples (divide by beta to get right number of 
% sample), PCSvec contains the probability of correct selection, in 
% expectation, given wvec, sval. 
%
% errstat is 0 if the value of sval is in a valid range for cfSoln, and is
% 1 if there was an error (meaning B() would be approximated by the
% approximate value function approximation in the cfSoln structure (e.g.,
% the KG* approximation)
%
% For purposes here, continuation region is approximated by the upper and
% lower boundary computed by the finite difference method, with the first
% order bias correction from Chernoff & Petkau methods, not the
% asymptotic approximations
%
% it is presumed that cfSoln was loaded by PDESolnLoad(), and that function
% in turn presumes that PDESolnCompute had been used to compute the PDE
% solution
%
% (c) 2015, 2016, S Chick, all rights reserved. Code is provided 'as is' with no
% warranty of correctness.
% UPDATED: 2021 04 19. Matlab seems to have changed behavior with interp2 -
% bombs now if checking a spot ok for s-coord but none in contin boundary
% in w-coord range - now does a check to keep crash from occuring.

%validate inputs, initialize some outputs
errstat = false;    % no error in computing B() (unless proven otherwise below
if (length(sval) > 1)
    warning('unpredictable behavior in PDEGetVals if sval is not a scalar');
    errstat = 1;
end
if nargin<4
    extraflag = false; % default: dont return PCS and E[N]
end
    
% Do a bit of data validation
% load in data from solution which will be required to compute the return
% values
stopfunc = cfSoln.Header.PDEparam.termrewardfunc;
approxfunc = cfSoln.Header.PDEparam.approxvaluefunc;
firsts = cfSoln.Header.firsts;
lasts = cfSoln.Header.lasts;
biggests = lasts(length(lasts));

wveccol=wvec(:); % insure we have a column vector for wvec available

% OUTPUT VALUE INITIALIZATION
% initialize Vvec, the output with value of continuing, to be the value of
% stopping immediately. 
stopvec = stopfunc(wveccol,sval,cfSoln.Header.PDEscale,cfSoln.Header.PDEparam);
Vvec = stopvec;
% initialize value of information of continuing to sample to 0, and PCS to
% be the PCS given stopping immediately... 
Bwsval = 0*Vvec;
if extraflag
    PCSvec = normpdf(abs(wveccol)/sqrt(sval),0,1);
% initialize the other outputs to the empty vector: their values will be
% filled in if they can be filled in below
    ENvec = 0*PCSvec;
% The values of the preceding outputs will be updated below if it is possible to do so
else
    PCSvec = [];
    ENvec = [];
end

% OUTPUT VALUE UPDATE IF POSSIBLE WITH SAVED FILES
% Find which block to use to estimate function value at time sval
if sval < firsts(1)  % sval was below the range for which data was originally computed
    % therefore, use approximation to value function provided during
    % computation time 
    Vvec = approxfunc(wveccol,sval,cfSoln.Header.PDEscale,cfSoln.Header.PDEparam);
    % from value function, subtract off reward from stopping immediately in
    % order to get the value of information from sequential sampling
    Bwsval = Vvec-stopvec; 
    % FIX: Might be able to do a smarter approximation by looking at
    % solution at smallest value of s which had computations.
    errstat = true;
%    warning( sprintf('sval = %f is below range of s, (%f, %f), computed: using approx func',sval,firsts(1),biggests));
elseif sval > biggests % sval is above the range for which data was originally computed
    Vvec = approxfunc(wveccol,sval,cfSoln.Header.PDEscale,cfSoln.Header.PDEparam);
    Bwsval = Vvec-stopvec;
    errstat = true;
%    warning( sprintf('sval = %f is above range of s, (%f, %f), computed: using approx func',sval,firsts(1),biggests));
else
    % find which block of the data structue we should look: put block
    % number into variable 'ind'
    isabove = (sval > lasts);
    [~, ind] = min(isabove);      % here, we know minimum will be 0 becase sval <= biggests, so b is index of file to use
    % get data for grid from block 'ind'
    wv = cfSoln.Data(ind).wvec;
    svec = cfSoln.Data(ind).svec;
    % interpolate upper and lower boundary from the bias-corrected boundary
    % from the finite difference.
    upper = interp1( svec, cfSoln.Data(ind).upvec, sval );
    lower = interp1( svec, cfSoln.Data(ind).downvec, sval );
    % interpolate the value of sampling at time sval for values of w in
    % wvec, for values in estimated continuation region
    if ~isempty(wveccol((wvec < upper) & (wvec > lower))) % SEC: 2021 04 19 FIX: bug fix for migration to 2021a - interp2 changed behavior when out of range for w vales on all matrix entries
        Bwsval((wvec < upper) & (wvec > lower))=interp2(svec,wv,cfSoln.Data(ind).Bwsmatrix,sval,wveccol((wvec < upper) & (wvec > lower)));
        %Bwsval((wvec < upper) & (wvec > lower))=interp2(svec,wv,cfSoln.Data(ind).Bwsmatrix,sval,wveccol((wvec < upper) & (wvec > lower)));,'linear',0);
        %Bwsval((wvec < upper) & (wvec >
        %lower))=interp2(svec,wv,cfSoln.Data(ind).Bwsmatrix,sval,wveccol((wvec
        %< upper) & (wvec > lower)));% SEC ORIG FROM PRE2020
    end
    if min(Bwsval) <= 0
        if min(Bwsval) < 0
            warning('interpolated value of continuing to sample computed to be negative');
        end
        VvecKG = approxfunc(wveccol,sval,cfSoln.Header.PDEscale,cfSoln.Header.PDEparam);
        % from value function, subtract off reward from stopping immediately in
        % order to get the value of information from sequential sampling
        BwsvalKG = VvecKG-stopvec; % get value of info from KG* lookahead
        Bwsval = max(Bwsval,BwsvalKG); % take max of PDE value and KG* type value
        errstat = true;
    end
    Vvec = Vvec + Bwsval;      % updated value for Vvec is value of continuing to sample plus value of stpping now
    if extraflag
        if ~isempty(wveccol((wvec < upper) & (wvec > lower))) % SEC: 2021 04 19 FIX: bug fix for migration to 2021a - interp2 changed behavior when out of range for w vales on all matrix entries
            ENvec((wvec < upper) & (wvec > lower)) = interp2(svec,wv,cfSoln.Data(ind).ENwsmatrix,sval,wveccol((wvec < upper) & (wvec > lower)));
        end %PCSvec interp2 does not use the check if in continuation region or not.
        PCSvec = interp2(svec,wv,cfSoln.Data(ind).EPCSwsmatrix,sval,wveccol);
    end
end


end

