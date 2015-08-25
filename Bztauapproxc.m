function [B]=Bztauapproxc(s, w0, s0);
% h is between 0 and s0
% s0 is a scalar for the reverse time brownian motion
% w0 = is the mean in the reverse time process
% April 2007 - "fixed" version
%tmp = sqrt( h / (1/tau0) / (1/tau0 + h));
%s 
%w0 
%s0
%1/s0 - 1/s
%exp(1/s0 - 1/s)
%sqrt(s0-s)
%PsiNorm( - w0 / sqrt(s0-s))
if (s > 0) * (s < s0)
%	if z0 > 0
%        B = exp(-h) * ( z0/tau0 + sqrt(tmp) * PsiNorm(- tmp * z0/tau0));
%        B = exp(-h) * ( z0/tau0 + tmp * PsiNorm( z0/tau0/tmp));
%	else
	B = exp(1/s0 - 1/s) * sqrt(s0-s) * PsiNorm( - w0 / sqrt(s0-s));
%	end
elseif s >= s0
%    'test'
    B = max(0,w0) - abs(s-s0);      % give dummy val for illegal s that is not as good
else
%    'test b'
    B = max(0,w0) - abs(s);      % give dummy val for illegal s that is not as good
end    
B = -B;      % Return - B (for use in minimization)
