% implements the psi function from Brezzi and Lai (2002), for use in
% approximating Gittins indices for multiarm bandit problems that are
% INFINITE HORIZON with DISCOUNTING
function [rvec]=bl_git_psi(svec,flag);
% svec is a vector of s = 1/tau values to check
% flag is true if a 'corrected' version is desired, false if the original version in Brezzi and Lai is desired.
[m n] = size(svec);
rvec=zeros(m,n);

if flag
	for i=1:m
        for j=1:n
            s=svec(i,j);
			if s < 0.1
                rval = sqrt(s/2);
			elseif s < 1
                rval = 0.445 - 0.07/sqrt(s);
			elseif s < 5
                rval = 0.62 - 0.245/sqrt(s);
			elseif s < 15
                rval = 0.7933 - 0.6327 /sqrt(s);
			else 
	            rval = sqrt(2 * log(s) - log (log (s)) - log (16 * pi) -.10573 ); % this one provides better fit
			end
            rvec(i,j)=rval;
        end
	end
else
	for i=1:m
        for j=1:n
            s=svec(i,j);
			if s < 0.2
                rval = sqrt(s/2);
			elseif s < 1
                rval = 0.49 - 0.11/sqrt(s);
			elseif s < 5
                rval = 0.63 - 0.26/sqrt(s);
			elseif s < 15
                rval = 0.77 - 0.58/sqrt(s);
			else 
                rval = sqrt(2 * log(s) - log (log (s)) - log (16 * pi)  );
	%            rval = sqrt(2 * log(s) - log (log (s)) - log (16 * pi) -.12 ); % this one provides better fit
			end
            rvec(i,j)=rval;
        end
	end
end
