
    % GET INDEX INTO RIGHT STUCTURE BASED ON TIME
    [a, b] = min(s0 > sbvec);
    ind=findex(b)
    wvec = cfSoln.Data(ind).wvec;
    svec = cfSoln.Data(ind).svec;
    Bwsmatrix = cfSoln.Data(ind).Bwsmatrix;
    Bmu0t0=mwig/scale.beta+interp2(svec,wvec,Bwsmatrix,s0,w0)/scale.beta
    Bw0s0=mwig + interp2(svec,wvec,Bwsmatrix,s0,w0)
    % This is the VOI info - need to add in the stopping to get the full
    % value function
end

PDEGetVals(cfSoln,wvec,sval) % try to get solution for values of w in wvec at time sval, given pde solution cfSoln



%%%%%%%%% convert to standardized version in (w,s) with w brownian motion
%%%%%%%%% in -s time scale
mwig = scale.beta * param.retire
w0= scale.beta*(mu0-param.retire)
s0 = 1 / (scale.gamma * param.t0)
sEND = 1 / (scale.gamma * param.tEND)

% set up info for plotting stuff
myfontsize=16
mysmallfontsize=14
points = 144*3 %spacing between labels on contours - made so that only one label appears per line


%[ta tb tc]=fminbnd('Bztauapproxc',0,1/scale.gamma/n0base,optimset('TolFun',scale.beta/scale.sigma,'MaxIter',10^4),0,1/scale.gamma/n0base);       % FIXED VERSION
[ta tb tc]=fminbnd('Bztauapproxc',0,s0,optimset('TolFun',scale.beta/scale.sigma,'MaxIter',10^4),0,s0);       % FIXED VERSION
-tb/scale.beta - scale.c/scale.gamma
sig0 = scale.sigma / sqrt(param.t0)
sig0*PsiNorm(-mu0/sig0) - (-tb/scale.beta - scale.c/scale.gamma)

%ds
%dt = 1/scale.gamma/s0 - 1/scale.gamma/(s0+ds)
%dw
%dmu=dw/scale.beta

%Bmu0t0=m+interp2(svec,wvec,Bwsmatrix,s0,w0)/scale.beta
%Bw0s0=scale.beta*m + interp2(svec,wvec,Bwsmatrix,s0,w0)
