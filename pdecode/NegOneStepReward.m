function [B]=NegOneStepReward(s, w0, s0, scale, param);
% s is between 0 and s0, so increment is from s0 down to s
% s0 is a scalar for the reverse time brownian motion
% w0 = is the mean in the reverse time process
% scale is the regular scaling, param gives the rest of the parameters
%1/s0 - 1/s
%exp(1/s0 - 1/s)
%sqrt(s0-s)
%PsiNorm( - w0 / sqrt(s0-s))

mu0 = w0/scale.beta;
t0 = 1 / (s0 * scale.gamma);
t = 1 / (s * scale.gamma);

if t > t0
    % first, figure out if this is a discounted reward or not
    if scale.discrate > 0 % we have discounting
        sampcosts = scale.c + (param.online==true)*mu0;
        sampcosts = sampcosts * (1 - scale.discrate^max(t-t0,0)) / (1 - scale.discrate);
        zvar = scale.sigma^2 * (t-t0) / (t0 * t);
    	B = exp(-scale.discrate*(t-t0)) * sqrt(zvar) * PsiNorm( - mu0 / sqrt(zvar));

        B2 = exp(1/s0 - 1/s) * sqrt(s0-s) * PsiNorm( - w0 / sqrt(s0-s));

    else % not discounted rewards
        sampcosts = (t-t0) * (scale.c + (param.online==true)*mu0);
        zvar = scale.sigma^2 * (t-t0) / (t0 * t);
    	B = -sampcosts + sqrt(zvar) * PsiNorm( - mu0 / sqrt(zvar));

        sampcosts2 = (1/s - 1/s0);
    	B2 = -sampcosts2 + sqrt(s0-s) * PsiNorm( - w0 / sqrt(s0-s));
    end
    
    sampcosts
    sampcosts2/scale.gamma
else
    B = mu0 - (t0 - t); % return B0 with zero samples, minus a penalty for t < t0
    B2 = w0 - (1/s - 1/s0)*scale.beta;
end
B
B2/scale.beta

B = -B;      % Return - B (for use in minimization)
