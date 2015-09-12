function [ pdescale, pdeparam ] = PDEInputConstructor( scalearray, paramarray  )
%PDEINPUTCONSTRUCTOR creates two data structures which are used as inputs
%to the free boundary solver code.
%
% Designed for computing solutions to free boundary problems of the sort
% found in Chick & Gans, Chick & Frazier, Chick Forster Pertile, Brezzi &
% Lai, etc. (variations on bayesian bandit problems).
%
% OUTPUTS: It returns two main parameters, one with 'scale' parameters,
% which has parameters from original problem definition (c, discrate,
% sigmaX, P=multiplier for adoption decision), and will hold values of
% scaled version of parameters for reverse time diffusion. 
%
% The other main output is paramarray, which has parameters which are used
% for the computation of outputs, file information, etc.
%
% The default fields for the arrays will be filled in. There is no need to
% call the Constructor with basicarray or advancedarry. However, if those
% parameters are passed, they are assumed to be a cell array of even
% length, with the names of the parameters and the default values assigned
% to it.
%       % basicarray = { 'c', 2.1, 'discrate', 0.0002 }, for example
%       It is ok to add extra parameters into the structures which are not
%       initially in the default outputs - these can be used in making
%       terminal reward functions more complex
%
% (c) 2015, S Chick
% Created: 2015 Aug 21
% Last touched: 2015 Aug 21
% 

PDELocalInit;

    % check if inputs were passed
    if nargin < 2
        paramarray = {};
    end
    if nargin < 1
        scalearray = {};
    end

    % set default values of parameters whose values can change without
    % needing to recalculate the PDE solution
    pdescale.c = 1;             % marginal cost per sample
    pdescale.sigma = 10e5;      % std deviation of sample
    pdescale.discrate = 0.0;    % continuous time discount rate per sample, e.g. 0.0 or .00001
    pdescale.P = 1;             % for inference of unknown mean W, reward on output is PW

    % set default values of parameters whose values can not be changed
    % without needing to recalculate the PDE solution
    pdeparam.termrewardfunc = @(wvec,s,p1,p2)max(wvec,0); % by default, reward is max(posterior mea, 0), will be scaled by pdescale.P
    %termrewardfunc: @(wvec,s,p1,p2)max(wvec,0)
    pdeparam.approxvaluefunc = @(wvec,s,p1,p2)max(wvec,0); % by default, reward is max(posterior mea, 0), will be scaled by pdescale.P
    %approxvaluefunc: @(wvec,s,p1,p2)PDECFApproxValue(wvec,s,p1)
    %approxvaluefunc: @(wvec,s,p1,p2)PDECGApproxValue(wvec,s,p1)
    pdeparam.approxmethod = @(s,p1,p2)CFApproxBoundW(s); % by default, reward is max(posterior mea, 0), will be scaled by pdescale.P
    %approxmethod: @(s,p1,p2)CGApproxBoundW(s)
    %approxmethod: [ .06 1500] % for CF, for example, vector gives the
    %initial value of s for the recursion, and the number of values of dw
    %to use in the grid
                                
    pdeparam.online = false;    % default: online learning, meaning results of patients in trial are counted in expected reward
    pdeparam.finiteT = false;   % default: false = infinite horizon problem, false = finite horizon problem
    pdeparam.t0 = 0.25;         % default: online learning, meaning results of patients in trial are counted in expected reward
    pdeparam.tEND = 20000;      % default: online learning, meaning results of patients in trial are counted in expected reward
    pdeparam.retire = 0;        % value of retirement option: if this is 0 or below it is ignored,
                                % retirement is also ignored when the
                                % discount rate is 0 (as in C&F), but is 
                                % used as the 'best alternative' if
                                % discount rate is positive (as in C&G)
    pdeparam.precfactor = 5;    % intended to be minimum number of grid points between 0 and upper bound of continuation set. should be min 2, bigger for finer grid
    pdeparam.DoPlot = true;     % true if diagnostic plots/text to be output during computations of pde
    pdeparam.figdir = PDEfigfilebase; %'Figure\';% directory name for figures to be output
    pdeparam.matdir = PDEmatfilebase; %'Matfiles\'; % directory name for storing .mat files with pde solutions

    % At present the following are declared but are not yet supported. not
    % needed for pde calculations, but might be useful for monte carlo at a
    % later time.
    pdeparam.UnkVariance = false;   % false if variance is known, true if unknown (default to known variance)
    pdeparam.DistributionType = []; % set to empty vector unless object used for sampling distribution
    pdeparam.Distribution = []; % 
%    pdescale.mu0 = 0;          % mean of prior distribution for uknown mean reward
%    pdeparam.UnkVariance = true;   % identify if the stopping boundary will be adjusted to account for sampling variance or not
%    pdeparam.DistributionType = @DistNormalMu; % identify the type of sampling distribution object
%    pdeparam.Distribution = advanced.DistributionType(pdescale.mu0, pdeparam.t0, pdescale.sigma); % create an instance of the object with the right hyperparameters
%
%    pdescale.mu0 = 0;          % mean of prior distribution for uknown mean reward
%    pdeparam.xi0 = 20;      % shape parameter for unknown variance
%    pdeparam.UnkVariance = true;   % identify if the stopping boundary will be adjusted to account for sampling variance or not
%    pdeparam.DistributionType = @DistNormalMuSig; % identify the type of sampling distribution object
%    pdeparam.Distribution = advanced.DistributionType(pdescale.mu0, pdeparam.t0, pdescale.sigma, pdeparam.xi0); % create an instance of the object with the right hyperparameters
%     pdeparam.UnkVarianceShape = -1.0; % Ignored for the moment. Might later be used for fudge factor for boundary

    pdeparam.fixedP = true; % set to true if expected reward on stopping is for P * expected reward per patient, false if patients not tested due to early stopping can also benefit from better alterantive
    
    [rval, pdescale2, pdeparam2] = PDEInputModifier(pdescale, pdeparam, scalearray, paramarray);

    if ~isfield(pdeparam2,'BaseFileName') % if no file name was passed, then come up with a default name based on discount rate value
        if pdescale2.discrate == 0.0  % continuous time discount rate per sample, e.g. 0.0 or .00001
            pdeparam2.BaseFileName = 'CF'; % base text for file names for output
% Should fix: get approxvaluefunc to take correct default if not finiteT
% and not passed in inputmodifier function
%            if ~pdeparam2.finiteT
%            end
    %approxvaluefunc: @(wvec,s,p1,p2)PDECFApproxValue(wvec,s,p1)
    %approxvaluefunc: @(wvec,s,p1,p2)PDECGApproxValue(wvec,s,p1)
        else
            pdeparam2.BaseFileName = 'CG'; % base text for file names for output
        end
    end
    pdeparam.approxvaluefunc = @(wvec,s,p1,p2)max(wvec,0); % by default, reward is max(posterior mea, 0), will be scaled by pdescale.P
    
    pdescale = pdescale2;
    pdeparam = pdeparam2;
% COMPUTED PARAMETERS: If more are put here, please fix the code to insure
% that they don't overwrite parameter values which may have been entered by
% the end user.
if ~rval
    pdescale = [];
    pdeparam = [];
end


end

