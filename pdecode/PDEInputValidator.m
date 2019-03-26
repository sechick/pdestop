function [ scaleout, paramout, rval, errorlist ] = PDEInputValidator( scale, param )
%PDEInputValidator checks the scale and param data structures for
%validity, and corrects values if needed, for usage in the Delay Differential 
% equations for the delay paper of Chick, Forster, Pertile (in alpha order). 
%It returns two main parameters, one with 'scale' parameters, and which 
%should be used instead of the scale and param parameters which are
%input to this routine. That is, use the scale and param OUT parameters
%for passing into DelayCurvesRecur (stage II of the backwards recursion, 
%which is a backwards recursion for a continuous time stopping problem),
%and then DelayStageOne (stage I of the backwards recursion, which is a
%one-stage lookahead story).
%
% INPUTS:
%   scale, param : two structures with parameters for the sequential
%       sampling problem and its computation. Typically this should be created
%       by calling the constructor method, PDEInputConstructor.
% OUTPUTS:
%   scaleout, paramout: updated copies of the structures which have been
%       passed as inputs, but with values 'tweaked' in order to make them
%       amenable for the computations. Some values are complicated
%       functions of the other parameters, and so it is important to use
%       the 'out' version of the parameters when calling the recursions.
%   rval: true if all is ok, false if there were some errors worth noting.
%   errorlist: text string which contains messages which explains some of
%       the errors worth noting, and how PDEInputValidator went about
%       fixing them.
%
% NOTE: This routine may create new fields
%
% Designed for Delay Sequential Trials paper of Chick, Forster, Pertile
% (alpha order).
%
% (c) 2015, S Chick, 
% Created: 2015 02 22
% Last touched: 2015 08 22
% 

routinename='PDEInputValidator';
errorlist = '';

[defaultscale, defaultparam] = PDEInputConstructor();  % find the default values used in the construtor method.

%currentDate = datestr(now,'mmmdd');
%myStruct.(currentDate) = [1,2,3]
rval = 1;

if scale.sigma <= 0
    scale.sigma = defaultscale.sigma;  % standard deviation of sample differences (between people in two trial arms)
    errorlist = sprintf('%s\n%s: sampled std dev was negative: sigma reset to %f\n',errorlist,routinename,scale.sigma);
    rval = 0;
end
if scale.c < 0
    scale.c = defaultscale.c;  % cost per sample
    errorlist = sprintf('%s\n%s: negative sampling cost, c reset to %f\n',errorlist,routinename,scale.c);
    rval = 0;
end
if (scale.discrate < 0.0)
    tmpval = scale.discrate;
    scale.theta = defaultscale.discrate;  % discount factor
    errorlist = sprintf('%s\n%s: invalid discount factor (%f), theta reset to %f\n',errorlist,routinename,tmpval,scale.discrate);
    rval = 0;
end
if scale.P < 0
    tmpval = scale.P;
    scale.P = defaultscale.P; % multiplier times unkonwn mean for payoff
    errorlist = sprintf('%s\n%s: Factor P should be positive, reset to %f\n',errorlist,routinename,tmpval);
    rval = 0;
end


%param.online = false;    % default: offline learning, meaning results of samples are not counted in expected reward
%param.finiteT = false;    % default: infinite horizon problem
if param.t0 < 0.001
    param.t0 = defaultparam.t0;  % effective number of samples in prior distribution for unknown mean
    errorlist = sprintf('%s\n%s: t0 should be at least 0.001, reset to %f\n',errorlist,routinename,param.t0);
    rval = 0;
end
if param.tEND < param.t0
    param.tEND = 2 * param.t0;  % end of time horizon should be after start of time horizon
    errorlist = sprintf('%s\n%s: tEND too small relative to t0, and was reset to %f\n',errorlist,routinename,param.tEND);
    rval = 0;
end
%param.retire = 0.0;    % value of retirement option (might be ignored in code for the moment)
if param.precfactor < 2
    param.precfactor = defaultparam.precfactor;  % end of time horizon should be after start of time horizon
    errorlist = sprintf('%s\n%s: precfactor should be at least 2, and was reset to %f\n',errorlist,routinename,param.precfactor);
    rval = 0;
end
if ~param.DoFileSave
    errorlist = sprintf('%s\n%s: WARNING: set DoFileSave to false at own risk\n',errorlist,routinename);
end    
%param.DoPlot = true;    % default: print diagnostic info during computation
if length(param.figdir) > 0
    if ~isdir(param.figdir)
        mkdir(param.figdir)
    end
end
if length(param.matdir) > 0
    if ~isdir(param.matdir)
        mkdir(param.matdir)
    end
end
if length(param.BaseFileName) == 0
    param.BaseFileName = defaultparam.BaseFileName;
    errorlist = sprintf('%s\n%s: base file name for pde solution output set to default, %f\n',errorlist,routinename,param.BaseFileName);
end

if ~isa(param.termrewardfunc, 'function_handle')
    param.termrewardfunc = defaultparam.termrewardfunc;
    errorlist = sprintf('%s\n%s: termrewardfunc should be a function handle, such as @(wvec,s,p1,p2)max(wvec,0).\n',errorlist,routinename);
end
if ~isa(param.approxvaluefunc, 'function_handle')
    param.approxvaluefunc = defaultparam.approxvaluefunc;
    errorlist = sprintf('%s\n%s: approxvaluefunc should be a function handle, such as @(wvec,s,p1,p2)max(wvec,0).\n',errorlist,routinename);
end
if ~isa(param.approxmethod, 'function_handle')
    if length(param.approxmethod) ~= 2
        param.approxmethod = defaultparam.approxmethod;
        errorlist = sprintf('%s\n%s: approxmeth should be a function handle to approximate upper stopping bound (e.g. CFApproxBoundW or CGApproxBoundW) or a vector of length 2 with [ dw numdw ] for starting value of s.\n',errorlist,routinename);
    end
end
% following are ignored for the moment
%param.UnkVariance
%param.DistributionType
%param.Distriubtion
%param.fixedP
%if length(param.DistributionType) > 0    % if a sampling distribution type is given, check to see if it is implemented in current directory, and check to see if the distribution instance is of the correct object type
% Following commented code identifies whether a given object is implemented
% or not in the 'current directory'.  removed, so that it does not require
% the sampling distribution object to be implemented in the current working
% directory - can be elsewhere on path
%    dirData = dir();      %# Get the data for the current directory
%    dirIndex = [dirData.isdir];  %# Find the index for directories
%    fileList = {dirData(dirIndex).name};  % Get list of directories
%    definedRVs={fileList{strncmp(fileList, '@Dist', 4)}};  % set of RVs available should include those in directories starting with '@Dist'
%    if ~sum(strcmp(definedRVs, strcat('@',func2str(param.DistributionType))))
%        errorlist = sprintf('%s\n%s: warning, param.DistributionType %s not found in current directory\n',errorlist,routinename,strcat('@',func2str(param.DistributionType)));        
%    end
%    if ~strcmp(class(param.Distribution),func2str(param.DistributionType))
%        errorlist = sprintf('%s\n%s: warning, param.DistributionType is %s but instance param.Distribution is of type %s\n',...
%            errorlist,routinename,strcat('@',func2str(param.DistributionType)), class(param.Distribution));        
%    end
%end

%rval = isempty(errorlist);
%if param.verbose && (~rval)
if ~isempty(errorlist)
    warning(errorlist)
end
if rval
    scaleout = PDEScaleStandardize(scale);
    paramout = param;
else
    scaleout = [];
    paramout = [];
end

end

