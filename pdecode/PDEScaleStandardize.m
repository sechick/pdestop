function PDEScaleOut = PDEScaleStandardize(PDEScale)
% PDEScaleStandardize: Takes in a structure with fields c, sigma, discrate,
% and P, and outputs alpha, beta, gamma, kappainv, according to
% transformations of C&F (when discrate = 0 and P=1), C&G (when discrate >
% 0 and P=1), or CFP (WSC, when P \neq 1).
% Factor is not handled in this routine: handle it separately by
% manipulating c and sigma and discrate elsewhere.

    [ST,I] = dbstack;   % get name of this function, in case it is needed.

    % get parameters for optimal stopping problem in (y,t) space
    c = PDEScale.c;             % cost per single sample
    sigma = PDEScale.sigma; 	% standard deviation of sample
    discrate = PDEScale.discrate;% (continuous time) discount rate per sample
    %P = PDEScale.P;             % multiplier times the mean reward for the adoption decision

    PDEScaleOut = PDEScale;

    if discrate == 0 % this is undiscounted case of C&F
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

end
