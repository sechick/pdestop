function PDEScaleOut = PDEScaleStandardize(PDEScale)
% PDEScaleStandardize: Takes in a structure with fields c, sigma, discrate,
% and P, and outputs alpha, beta, gamma, kappainv, according to
% transformations of C&F (when discrate = 0 and P=1), C&G (when discrate >
% 0 and P=1), or CFP (WSC, when P \neq 1).

[ST,I] = dbstack;   % get name of this function, in case it is needed.

% get parameters for optimal stopping problem in (y,t) space
c = PDEScale.c;             % cost per single sample
sigma = PDEScale.sigma; 	% standard deviation of sample
discrate = PDEScale.discrate;% (continuous time) discount rate per sample
P = PDEScale.P;             % multiplier times the mean reward for the adoption decision

PDEScaleOut = PDEScale;

if P == 1
    if discrate == 0 %
        PDEScaleOut.alpha = c^(1/3) * sigma^(-4/3);
        PDEScaleOut.beta = c^(-1/3) * sigma^(-2/3);
        PDEScaleOut.gamma = c^(2/3) * sigma^(-2/3);
        PDEScaleOut.kappainv = 0;                   % this value is unused when disc rate is 0   
    else % set up parameters for discounted case of C&G
        PDEScaleOut.alpha = discrate^(1/2) * sigma^(-1);
        PDEScaleOut.beta = discrate^(-1/2) * sigma^(-1);
        PDEScaleOut.gamma = discrate;
        PDEScaleOut.kappainv = c * discrate^(-3/2) * sigma^(-1);    
    end
else
    % FIX: need to look at standardized problem and make sure that P<>1 is
    % handled both here and delayoffinesimplereward both. Need also to
    % handle case of I > 0 externally by tweaking the means, but that is
    % separate from the standardization.
    if discrate == 0 %
        c = c/P;
        PDEScaleOut.alpha = c^(1/3) * sigma^(-4/3);
        PDEScaleOut.beta = c^(-1/3) * sigma^(-2/3);
        PDEScaleOut.gamma = c^(2/3) * sigma^(-2/3);
        PDEScaleOut.kappainv = 0;                   % this value is unused when disc rate is 0   
    else % set up parameters for discounted case of C&G
        sigma = sigma/P;
%        c = c/P;
        PDEScaleOut.alpha = discrate^(1/2) * sigma^(-1);
        PDEScaleOut.beta = discrate^(-1/2) * sigma^(-1);
        PDEScaleOut.gamma = discrate;
        PDEScaleOut.kappainv = c * discrate^(-3/2) * sigma^(-1);    
    end
    warning('FIX: need to confirm %s for P <> 1 and both c, discrate > 0',ST.name);
end

end
