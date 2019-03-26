function [rval, pdescaleout,pdeparamout] = PDEInputModifier(pdescale, pdeparam, scalearray, paramarray)
%PDEInputModifier creates two data structures which are used as inputs
%to the free boundary solver code.
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
% call the Constructor with scalearray or paramarray. However, if those
% parameters are passed, they are assumed to be a cell array of even
% length, with the names of the parameters and the default values assigned
% to it.
%       % scalearray = { 'c', 2.1, 'discrate', 0.0002 }, for example
%
% (c) 2015, S Chick
% Created: 2015 Aug 21
% Last touched: 2015 Aug 21

if nargin < 4
    paramarray = {};
end
if nargin < 3
    scalearray = {};
end

scalearrayarrlen = length(scalearray);
rval = 1;
for i=1:(scalearrayarrlen/2)
    if isfield(pdescale,scalearray{2*i-1}) %|| strcmp(scalearray{2*i-1},'mumax') ||  strcmp(scalearray{2*i-1},'mumin')
        pdescale.(scalearray{2*i-1}) = scalearray{2*i};
    elseif isfield(pdeparam,scalearray{2*i-1})
        warning('attempt ignored to add field of param to scale structure: %s',char(scalearray{2*i-1})); 
%        pdescale.(scalearray{2*i-1}) = scalearray{2*i};
        rval = 0;
    else
        fprintf('adding pdescale field: %s, with value %s\n',char(scalearray{2*i-1}),char(scalearray{2*i})); 
        pdescale.(scalearray{2*i-1}) = scalearray{2*i};
%        rval = 0;
    end
end

paramarrlen = length(paramarray);
for i=1:(paramarrlen/2)
    if isfield(pdeparam,paramarray{2*i-1}) % || ...
%            strcmp(paramarray{2*i-1},'dmu') || ...
 %           strcmp(paramarray{2*i-1},'simFreqDeltaVec') 
        pdeparam.(paramarray{2*i-1}) = paramarray{2*i};
    elseif isfield(pdescale,paramarray{2*i-1})
        warning('attempt ignored to add field of scale to param structure: %s',char(scalearray{2*i-1})); 
%        pdeparam.(paramarray{2*i-1}) = paramarray{2*i};
        rval = 0;
    else
        fprintf('adding pdeparam field: %s with value %s\n',char(paramarray{2*i-1}),char(paramarray{2*i})); 
        pdeparam.(paramarray{2*i-1}) = paramarray{2*i};
%        rval = 0;
    end
end

pdescaleout = pdescale;
pdeparamout = pdeparam;

end

