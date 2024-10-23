function x = randmn(p,n)
% Generate samples from a multinomial distribution.
% Input:
%   p: k dimensional probability vector
%   n: number of samples
% Ouput:
%   x: 1 by n generated samples x~Mul(p)
if nargin == 1
    n = 1;
end

r = rand(1,n);
p = cumsum(p(:));
[~,x] = histc(r,[0;p/p(end)]);

return;