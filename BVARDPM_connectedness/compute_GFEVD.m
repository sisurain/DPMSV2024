function [GFEVD, Theta_unnorm] = compute_GFEVD(B_lags, Sigma, H)
% =========================================================================
% compute_GFEVD  Generalized Forecast Error Variance Decomposition
%
% Computes the H-step generalized forecast error variance decomposition
% of Pesaran & Shin (1998), with row normalization following Diebold &
% Yilmaz (2012, IJF) and Chan & Yu (2022, ICEEE) Section 6.1.
%
% MODEL
%   y_t = B_1 y_{t-1} + ... + B_p y_{t-p} + u_t,    u_t ~ N(0, Sigma)
%
%   The function takes the REDUCED-FORM VAR coefficients B_1, ..., B_p
%   and the (instantaneous) reduced-form covariance Sigma. Whether Sigma
%   is time-varying or constant is the caller's responsibility:
%
%   - For DPM-OISV: Sigma_t = lambda_t * B0^{-1} * D_t * (B0^{-1})'
%                   computed inside the MCMC loop, passed in here as Sigma
%   - For CCMM SV : Sigma_t = A^{-1} * Lambda_t * (A^{-1})'
%                   same idea, just with the unit-lower-triangular A
%
%   This function is agnostic to the source of Sigma. It only requires that
%   Sigma is symmetric positive-definite.
%
% PESARAN-SHIN GENERALIZED FEVD
%
%   The unnormalized contribution of variable j to the H-step ahead forecast
%   error variance of variable i is:
%
%       theta_{ij}^g(H) = sigma_{jj}^{-1} * sum_{h=0}^{H-1} (e_i' A_h Sigma e_j)^2
%                       / sum_{h=0}^{H-1} (e_i' A_h Sigma A_h' e_i)
%
%   where A_h are the reduced-form moving-average coefficients (A_0 = I_n),
%   computed by iterating the VAR forward:
%       A_0 = I_n
%       A_h = sum_{l=1}^{min(h,p)} B_l * A_{h-l}     for h >= 1
%
%   Following DY (2012) and Chan-Yu (2022), the rows are then normalized:
%       C_{i<-j}(H) = theta_{ij}^g(H) / sum_k theta_{ik}^g(H)
%
%   so each row sums to 1. C_{i<-j} measures the share of variable i's
%   H-step forecast error variance attributable to shocks to variable j.
%
% INPUTS
%   B_lags : n x n x p   array of reduced-form VAR coefficient matrices
%                        B_lags(:,:,l) is B_l (the lag-l coefficient).
%                        IMPORTANT: these are the REDUCED-FORM coefficients,
%                        i.e., the regression matrix mapping y_{t-l} -> y_t,
%                        NOT structural-form coefficients.
%   Sigma  : n x n       symmetric positive-definite reduced-form covariance
%   H      : scalar      forecast horizon (>= 1). Standard DY value is 10.
%
% OUTPUTS
%   GFEVD         : n x n   row-normalized GFEVD matrix.
%                           GFEVD(i,j) = share of variable i's H-step FEV
%                           attributable to shocks to variable j.
%                           Each row sums to 1.
%   Theta_unnorm  : n x n   unnormalized GFEVD (before row normalization).
%                           Useful for diagnostics; row sums need not equal 1.
%
% USAGE EXAMPLE
%   B = ...; Sigma = ...;     % from one MCMC draw
%   GFEVD = compute_GFEVD(B, Sigma, 10);
%
%   % Total connectedness (Diebold-Yilmaz system-wide measure):
%   C_total = (sum(GFEVD(:)) - trace(GFEVD)) / n * 100;
%
%   % "From others" for variable i (row sum minus diagonal):
%   from_i = (sum(GFEVD(i,:)) - GFEVD(i,i)) / n * 100;
%
%   % "To others" for variable j (column sum minus diagonal):
%   to_j = (sum(GFEVD(:,j)) - GFEVD(j,j)) / n * 100;
%
% REFERENCES
%   Pesaran, H. and Shin, Y. (1998). "Generalized impulse response analysis
%       in linear multivariate models." Economics Letters 58, 17-29.
%   Diebold, F. and Yilmaz, K. (2012). "Better to give than to receive:
%       Predictive directional measurement of volatility spillovers."
%       International Journal of Forecasting 28, 57-66.
%   Diebold, F. and Yilmaz, K. (2014). "On the network topology of variance
%       decompositions: Measuring the connectedness of financial firms."
%       Journal of Econometrics 182, 119-134.
%   Chan, J. and Yu, X. (2022). "Fast and accurate variational inference
%       for large Bayesian VARs with stochastic volatility." Section 6.1.
%
% =========================================================================

%% ----- Input validation -----
if nargin < 3
    error('compute_GFEVD:NotEnoughInputs', ...
        'Required: B_lags (n x n x p), Sigma (n x n), H (scalar).');
end

[n1, n2, p] = size(B_lags);
if n1 ~= n2
    error('compute_GFEVD:NonSquareB', ...
        'Each lag matrix in B_lags must be square; got %d x %d.', n1, n2);
end
n = n1;

if any(size(Sigma) ~= [n n])
    error('compute_GFEVD:SigmaDimMismatch', ...
        'Sigma must be %d x %d to match B_lags; got %d x %d.', ...
        n, n, size(Sigma,1), size(Sigma,2));
end

if ~isscalar(H) || H < 1 || H ~= round(H)
    error('compute_GFEVD:BadHorizon', ...
        'H must be a positive integer scalar; got %s.', mat2str(H));
end

% Symmetrize Sigma (defensive against floating-point asymmetry)
Sigma = 0.5 * (Sigma + Sigma');

%% ----- Compute reduced-form MA coefficients A_0, A_1, ..., A_{H-1} -----
%
% Recursion: A_0 = I, A_h = sum_{l=1}^{min(h,p)} B_l * A_{h-l}
%
% Stored as a 3-D array A_h(:,:,h+1) for h = 0, 1, ..., H-1
%   (MATLAB indexing: A_h(:,:,1) is A_0, A_h(:,:,2) is A_1, etc.)
A_h = zeros(n, n, H);
A_h(:,:,1) = eye(n);
for h = 1 : H-1
    Ah = zeros(n, n);
    for l = 1 : min(h, p)
        Ah = Ah + B_lags(:,:,l) * A_h(:,:,h-l+1);
    end
    A_h(:,:,h+1) = Ah;
end

%% ----- Compute Theta_unnorm (Pesaran-Shin generalized variance shares) -----
%
% Numerator   : sum_{h=0}^{H-1} (e_i' A_h Sigma e_j)^2  / sigma_{jj}
% Denominator : sum_{h=0}^{H-1} (e_i' A_h Sigma A_h' e_i)
%
% We vectorize this entire computation by accumulating the n x n outer
% products M_h := A_h * Sigma over h.
%
% Then:
%   numerator(i,j)   = (1/sigma_{jj}) * sum_h M_h(i,j)^2
%   denominator(i)   = sum_h (M_h * A_h')(i,i)
%                    = sum_h sum_k A_h(i,k) * Sigma(k,:) * A_h(i,:)'
%
% The denominator is variable-i-only (does not depend on j), so we compute
% it as an n-vector and then broadcast across columns.

sigma_jj = diag(Sigma);             % n x 1: diagonal of Sigma

% Numerator: sum_h (M_h .^ 2), where M_h = A_h * Sigma
num_acc = zeros(n, n);
% Denominator: sum_h diag(M_h * A_h') = sum_h diag(A_h * Sigma * A_h')
den_acc = zeros(n, 1);

for h = 1 : H
    Ah = A_h(:,:,h);                % n x n
    Mh = Ah * Sigma;                % n x n: Mh(i,j) = e_i' A_h Sigma e_j

    num_acc = num_acc + Mh.^2;

    % diag(Mh * Ah') without forming the full n x n product:
    %   (Mh * Ah')(i,i) = sum_k Mh(i,k) * Ah(i,k)
    %                   = sum across columns of element-wise product
    den_acc = den_acc + sum(Mh .* Ah, 2);
end

% Numerator after sigma_jj normalization: divide each column j by sigma_jj
%   Theta_unnorm(i,j) = num_acc(i,j) / sigma_jj(j) / den_acc(i)
% Use bsxfun-style broadcasting:
Theta_unnorm = num_acc ./ (den_acc * sigma_jj');   % n x n

%% ----- Row normalization (Diebold-Yilmaz convention) -----
% Each row of GFEVD sums to 1 by construction.
row_sums = sum(Theta_unnorm, 2);            % n x 1
% Defensive: replace any zero row-sums (which should not occur for PD Sigma)
% to avoid divide-by-zero. If a zero row-sum does occur, set that row to
% 1/n uniformly so the rest of the matrix is still usable.
zero_rows = (row_sums <= 0) | ~isfinite(row_sums);
if any(zero_rows)
    warning('compute_GFEVD:ZeroRowSum', ...
        '%d row(s) have non-positive Theta_unnorm row sum; setting to 1/n.', ...
        sum(zero_rows));
    row_sums(zero_rows) = 1;
end

GFEVD = Theta_unnorm ./ row_sums;

% Force exact row-stochastic (defensive against floating-point drift)
% Only apply to non-degenerate rows.
GFEVD(zero_rows, :) = 1 / n;

end
