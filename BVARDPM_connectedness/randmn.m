function x = randmn(p, n)
% Generate samples from a categorical distribution using Inverse Transform Sampling.
% Input:
%   p: k dimensional probability vector (or weights)
%   n: number of samples (default 1)
% Output:
%   x: 1 by n generated samples (indices of the categories)

if nargin == 1
    n = 1;
end

% Generate uniform random numbers
r = rand(1, n);

% Calculate cumulative probabilities and ensure p is a column vector
p_cumsum = cumsum(p(:));

% Define the edges for discretization. Normalize and prepend 0.
edges = [0; p_cumsum / p_cumsum(end)];

% Robustness: Ensure the last edge is slightly above 1 to robustly capture 
% cases where rand() returns exactly 1. eps(1) is the smallest distance 
% from 1.0 to the next double-precision number.
edges(end) = 1 + eps(1);

% Use discretize to find the bin index for each random number r.
% The output x will have the same shape as r (1 by n).
x = discretize(r, edges);
end