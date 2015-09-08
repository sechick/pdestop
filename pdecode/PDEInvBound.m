function [ svec ] = PDEInvBound( bndfunchandle, wvec )
%PDEInvBound: invert the boundary in bndfunchande at the values in wvec, returning svec = bndfunchandle^{-1}(wvec).
%   Used for inverting the upper stopping boundary described by the function handle.

    svec = zeros(size(wvec)); % preallocate memory, by default value is 0 (to handle cases when svec(i) is incorrectly passed as <0)

    sqrdiff = @(x,y) (x-y)^2;

    lasts=1; % initial search point for inverse
    for i=1:length(wvec)
        if wvec(i) > 0
            svec(i) = fminsearch(@(x) sqrdiff(bndfunchandle(x), wvec(i)), lasts);
        end
        lasts = svec(i); % work on assumption that inverse of next point is near to inverse of current point, to hopefully speed search
    end

end

