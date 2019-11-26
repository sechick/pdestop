function [ evivec, numvec ] = DelayOfflineSimpleReward( wvec, s0, scale, param )
%DelayOfflineSimpleReward:  find standardized reward in (w, s) coordinates
%when there is a delayed observation .

    if isfield(param,'tau')
        tau = param.tau;
    elseif isfield(scale,'tau')
        tau = scale.tau;
    else
        tau = 0;
    end

if 0 % commented out: we need to handle P in the 'pdescalestandardize' and the 'pdegetvals' routinesmv De
    if isfield(param,'P')
        P = param.P;
    elseif isfield(scale,'P')
        P = scale.P;
    else
        P = 1;
    end
end

    timet = 1 / (scale.gamma * s0);                 % get time t equivalent of the time to test, s0.
    uinpipe = max(0, min(tau, timet - param.t0));   % compute number of patients in pipeline. 
    suinpipe = 1 / (scale.gamma * (timet + uinpipe)); % compute time when decision is made after pipeline patients seen, 
    
    sincrem = abs( suinpipe - s0);      % want to compute time accounting for fact num pipeline might be less than tau if in 'stage I' sampling
%    sincrem = s0 * (tau * scale.gamma * s0) / (1 + tau * scale.gamma * s0); % find ds in going from t=1/(gamma s0) to t+tau, in s coords

    %dincrem = sqrt(3*sincrem/2);    % implied reward from diffusion over that increment
    evivec = max(0, wvec); 
    if scale.discrate == 0
        evitmp = -(1/(s0-sincrem) -1/s0) + sqrt(sincrem) * PsiNorm(-wvec/sqrt(sincrem));
    else
        evitmp = exp(-(1/(s0-sincrem) -1/s0)) * sqrt(sincrem) * PsiNorm(-wvec/sqrt(sincrem));
    end
    evivec = max(evivec, evitmp);
    numvec = 0 * evivec;
end