function [ reward, numsamps ] = PDEsimplereward( wvec )
%PDEsimplereward: returns expected reward and 0 added samples to achieve it,
%   using the (w, s) time scale. Useful as generic terminal condition.

    reward = max(wvec, 0);
    numsamps = 0*reward;

end

