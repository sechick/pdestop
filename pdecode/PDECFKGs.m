function [ kgfactor, dsvec ] = PDECFKGs( wvec, s0, scale )
%PDECFKGs: Find the KGstar allocation, putting the kgfactor-1 in the output
%vector kgfactor, and the increment in s coordinates which gives that
%factor, in the standardized (w,s) reverse time scale, for the stopping
%problem with ZERO discount rate. This generalizes the KG* of Frazier and
%Powell, with extension to do a search on a grid of s values between the
%lower and upper bound given in their 2010 Decision Analysis paper, and
%with the generalization that it is the expected voi per unit cost of
%computation (rather than per unit sample), as we allow for samples from
%different alternatives to have different costs, scale.c.
%
%In this (w,s) scaling, it is important to know the conversion factor gamma
%(in the language of C&F) for the reverse time diffusion, as this
%determines the increment in s associated with taking 1 sample in the
%original problem: this is to be passed as a field in 'scale' or as the
%value of scale itself.
%
% If scale is not passed as an argument, then a default vector of values of
% svec is created. If scale is passed as a real, then it is interpreted as
% gamma, the rescacing parameter s = 1/(gamma n_t), for the -s scale
% diffusion. If scale is passed as a structure, then the .gamma field is
% used for the value of gamma.
%
% The value of the optimal one-step lookahead, in -s coordinates, is returned in
% dsvec. This would convert to a lookahead of 1/(gamma ds) samples in t
% coordinates
%
% At worst case, assume that one can stop immediately and get the best of 0
% and the value of wvec. This is the usual terminal reward for stopping.

kgfactor = - ones(size(wvec)); % initialize the kgfactor to be a low value
dsvec = 0*kgfactor;    % initialize the increment for the optimal sampling plan to be 0: it will increase above 0 in the continuation region in calculations below.

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
        kgtmp = - 1 + sqrt(sincrem(i)) * PsiNorm(-abs(wvec)/sqrt(sincrem(i))) / (1/svec(i) - 1/s0);
        dsvec(kgtmp > kgfactor) = sincrem(i);
        kgfactor = max(kgfactor,kgtmp);
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
    svec = s0 * 1 ./ ( 1 + tstincrem*gamma*s0)
    sincrem = s0-svec(svec>0 & svec < s0) ;            % time elapse from scur to valid values in svec    
%    figure(1001); hold off; figure(1002); hold off; figure(1003); hold off;
    for i=1:length(sincrem)           % check the lookaheads over that interval
        kgtmp = - 1 + sqrt(sincrem(i)) .* PsiNorm(-abs(wvec)/sqrt(sincrem(i))) ./ (1./(s0-sincrem(i)) - 1/s0);
        dsvec(kgtmp > kgfactor) = sincrem(i) ;
%        costtmp = -(1/(s0-sincrem(i)) - 1/s0);
%        costvec(evitmp > evivec) = costtmp;
        kgfactor = max(kgfactor,kgtmp);
        figure(1001); plot(wvec,sincrem(i)); hold on; plot(wvec,s0);
        figure(1002); plot(wvec,evivec,'k',wvec,evitmp,'b'); hold on;
        figure(1003); plot(wvec,costvec); hold on;
    end
    
    % now check bounds from Frazier and Powell 2010
    kgfactor2 = max(wvec,0);   % lower bound on best you can do is to stop immediately
    dsvec2 = 0*kgfactor2;   % lower bound on best you can do is to stop immediately
    NUMKGCHECKS = 12;
    CHECKVEC=6:(4/NUMKGCHECKS):10;
    kgstarterm = wvec.^2 / s0; 
    figure(10001); hold off; figure(10002); hold off;
    for i=1:length(CHECKVEC)
        sfactor = ( 1 /(4*gamma*s0)  ) * (kgstarterm - 1 + sqrt(kgstarterm.^2 + CHECKVEC(i)*kgstarterm + 1)); % m-underbar from prop 4, frazier & powell, 2010
        sincrembase = s0 * max( gamma*s0 / (1 + gamma*s0) , sfactor*gamma*s0 ./ ( 1 + sfactor*gamma*s0) ); % this is the increment for kgstar for C&F. 
        kgtmp = - 1 + sqrt(sincrembase) .* PsiNorm(-abs(wvec)./sqrt(sincrembase)) . (1./(s0-sincrembase) - 1/s0);
        dsvec2(kgtmp > kgfactor2) = sincrembase(kgtmp > kgfactor2);%(evitmp > evivec2);
        kgfactor2 = max(kgfactor2,kgtmp);
        figure(10001); plot(wvec,sincrembase); hold on; plot(wvec,s0);
        figure(10002); plot(wvec,kgfactor2,'k',wvec,evitmp,'b'); hold on;
%        figure(10003); plot(wvec,costvec); hold on;    end
    kgfactor-kgfactor2;
    figure(2000)
    plot(wvec,evivec,'-',wvec,kgfactor2,'--')
    figure(2001)
    plot(wvec,dsvec,'-',wvec,dsvec2,'--')

    kgfactor = max(kgfactor,kgfactor2);
    dsvec = max(dsvec,dsvec2);

end


end
