% sample_SV_OISV.m
function [h, S] = sample_SV_OISV(ystar, h, rho, sig2, mu)
% Samples the log-volatilities using the KSC (1998) algorithm.
% Optimized for efficiency (vectorization) and numerical robustness.
% Inputs:
%   ystar: T x 1 vector of log-squared standardized residuals (log(s2))
%   h: T x 1 vector of current log-volatilities
%   rho: Scalar persistence parameter (phi in main script)
%   sig2: Scalar innovation variance (omega2 in main script)
%   mu: Scalar unconditional mean (0 in OISV)

T = length(h);
n_mix = 7;

% Ensure inputs are column vectors for consistency
if isrow(ystar); ystar = ystar'; end
if isrow(h); h = h'; end

% 1. KSC 7-component normal mixture parameters (Persistent Constants)
persistent m_N sig2_N log_p_N log_const_N
if isempty(m_N)
    p_N = [0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750];
    % Means adjusted by E[log(chi^2(1))] = -1.2704
    m_N = [-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819] - 1.2704;
    sig2_N = [5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261];
    % Pre-calculate constants for the log-pdf calculation
    log_p_N = log(p_N);
    % log_const_N = -0.5 * log(2*pi*sig2_N)
    log_const_N = -0.5 * (log(2*pi) + log(sig2_N));
end

% 2. Sample Mixture Components S (Vectorized and Robust Log-Scale)

% Calculate the mean of the observation equation for each component: h + m_N
% Use broadcasting (implicit expansion) instead of repmat
means_mix = h + m_N; % (T x 7)

% Calculate log-likelihood (log_L) for each observation under each component
% log_L = log_const_N - 0.5 * (ystar - mean)^2 / sig2_N
log_L = log_const_N - 0.5 * ((ystar - means_mix).^2 ./ sig2_N); % (T x 7)

% Calculate log posterior probabilities: log(q) = log(p_N) + log(L)
log_q = log_p_N + log_L;

% Normalize probabilities using Log-Sum-Exp trick for stability
max_log_q = max(log_q, [], 2);
% Subtract max before exponentiating to prevent overflow/underflow
q = exp(log_q - max_log_q);
q = q ./ sum(q, 2); % Normalize rows to sum to 1

% Sample S from the discrete distribution q
tmprand = rand(T, 1);
% Efficient sampling using cumulative sum comparison
S = n_mix - sum(tmprand < cumsum(q, 2), 2) + 1; % (T x 1)

% 3. Sample h using the Precision Sampler

% Extract parameters based on the sampled states S
d = m_N(S)'; % (T x 1) vector of offsets
iOmega_vec = 1./sig2_N(S)'; % (T x 1) vector of inverse observation variances

% Construct the State Precision Matrix (Q = HiSH/sig2)
% This matrix is derived from the stationary AR(1) process structure.
% We construct it efficiently using spdiags.

iSig2 = 1/sig2;

% Main diagonal: [1, 1+rho^2, ..., 1+rho^2, 1] / sig2
diag_main = ones(T, 1);
if T > 1
    % For T >= 2. If T=2, (2:T-1) is empty, resulting in [1; 1].
    diag_main(2:T-1) = 1 + rho^2;
end
% Off-diagonal: -rho / sig2
diag_off = -rho * ones(T, 1);

% Construct sparse matrix Q
Q = spdiags([diag_off, diag_main, diag_off], -1:1, T, T) * iSig2;

% Full Precision Matrix Kh = Q + iOmega
% We add iOmega_vec directly to the diagonal of Q
Kh = Q + spdiags(iOmega_vec, 0, T, T);

% Calculate Posterior Mean h_hat
% Kh * h_hat = Q * (mu*1_T) + iOmega * (ystar - d)

% Contribution from data (using element-wise multiplication for the diagonal iOmega)
data_contrib = iOmega_vec .* (ystar - d);

% Contribution from prior mean mu
if abs(mu) > 1e-8
    % If mu != 0
    prior_mean_contrib = Q * (mu*ones(T,1));
    h_hat = Kh \ (prior_mean_contrib + data_contrib);
else
    % If mu=0 (OISV case), the calculation simplifies
    h_hat = Kh \ data_contrib;
end

% Draw h
% Robustness: Attempt Cholesky decomposition.
[L, flag] = chol(Kh, 'lower');

if flag == 0
    % Success: Draw using the Cholesky factor of the precision matrix
    % h = h_hat + inv(L')*Z
    h = h_hat + L' \ randn(T, 1);
else
    % Failure: Kh is not numerically positive definite.
    % Throw an error; the calling MCMC script catches it and keeps
    % the previous draw of h.
    error('sample_SV_OISV:CholeskyFailed', 'Posterior precision matrix Kh is not positive definite.');
end
end