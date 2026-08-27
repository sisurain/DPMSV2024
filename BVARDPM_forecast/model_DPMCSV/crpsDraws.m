function crps = crpsDraws(y, Xdraws)
% CRPSDRAWS  CRPS from predictive draws (Krueger et al. 2021)
%   crps = crpsDraws(y, Xdraws)
%   y      : scalar realized value
%   Xdraws : vector of draws from the predictive distribution

Xdraws = Xdraws(:);
Ndraws = length(Xdraws);

if any(isnan(Xdraws)) || isnan(y)
    crps = NaN; return;
end

Xordered = sort(Xdraws);
[alpha, beta] = deal(NaN(Ndraws-1,1));

ndx1 = y < Xordered(1:end-1);
ndx2 = (Xordered(1:end-1) <= y) & (y <= Xordered(2:end));
ndx3 = Xordered(2:end) < y;

alpha(ndx1) = 0;
foo = y - Xordered(1:end-1); alpha(ndx2) = foo(ndx2);
foo = Xordered(2:end) - Xordered(1:end-1); alpha(ndx3) = foo(ndx3);

foo = Xordered(2:end) - Xordered(1:end-1); beta(ndx1) = foo(ndx1);
foo = Xordered(2:end) - y; beta(ndx2) = foo(ndx2);
beta(ndx3) = 0;

ALPHA = zeros(Ndraws+1,1);
BETA  = ALPHA;
ALPHA(2:end-1) = alpha;
BETA(2:end-1)  = beta;
if y >= Xordered(end), ALPHA([1 end]) = y - Xordered(end); end
if y <= Xordered(1),   BETA([1 end])  = Xordered(1) - y;   end

ndx  = (0:Ndraws)' / Ndraws;
crps = sum(ALPHA .* ndx.^2 + BETA .* (1 - ndx).^2);
end
