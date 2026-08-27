% sample_SV0para_OISV.m
function [phi, sig2, flag_phi] = sample_SV0para_OISV(h, phi, Hyper)
% Samples the SV parameters phi (persistence) and sig2 (innovation variance)
% for a zero-mean AR(1) process. Fully vectorized and numerically robust.
%
% Inputs:
%   h: T x n matrix of log-volatilities
%   phi: n x 1 vector of current persistence parameters
%   Hyper: Structure containing priors (nuh, Sh, phi0, Vphi - all n x 1)

[T, n] = size(h);

% Ensure phi is a column vector for consistency
if isrow(phi); phi = phi'; end
h1_sq = h(1,:)'.^2; % Squared initial states (n x 1)

% Define stability constraints
STATIONARITY_BOUND = 0.999;
MIN_SIG2 = 1e-12;

% 1. Sample sig2 (Gibbs Step) - Vectorized

% Calculate innovations e_h efficiently
% Innovations from t=2 to T (using broadcasting)
e_h_innov = h(2:end,:) - h(1:end-1,:) .* phi'; % (T-1) x n
% Contribution of initial state h1 (scaled by sqrt(1-phi^2))
% SSE = sum(e_h_innov^2) + h1^2 * (1-phi^2)
SSE = sum(e_h_innov.^2, 1)' + h1_sq .* (1 - phi.^2); % (n x 1)

% Posterior parameters (Inverse Gamma)
shape_post = Hyper.nuh + T/2;
scale_post = Hyper.Sh + SSE/2;

% Sample sig2. IG(A, B) -> 1/Gamma(A, 1/B)
sig2 = 1./gamrnd(shape_post, 1./scale_post); % (n x 1)

% Robustness: Enforce minimum variance
sig2(sig2 < MIN_SIG2) = MIN_SIG2;
iSig2 = 1./sig2; % Inverse variance

% 2. Sample phi (Vectorized Metropolis-Hastings Step)

% Calculate proposal parameters (Gaussian approximation)
h_lag = h(1:T-1, :);
h_lead = h(2:T, :);

% Kphi (Precision) = 1/Vphi + sum(h_lag^2)/sig2
Kphi = 1./Hyper.Vphi + sum(h_lag.^2, 1)' .* iSig2;
% phi_hat (Mean) = (phi0/Vphi + sum(h_lag*h_lead)/sig2) / Kphi
phi_hat = (Hyper.phi0./Hyper.Vphi + sum(h_lag .* h_lead, 1)' .* iSig2) ./ Kphi;

% Draw candidate phic
phic = phi_hat + (1./sqrt(Kphi)) .* randn(n, 1);

% 3. Calculate Acceptance Probability (Vectorized and Robust)

% The MH ratio depends on the adjustment term g(phi) due to the initial condition:
% g(phi) = 0.5*log(1-phi^2) - 0.5*h1^2*(1-phi^2)/sig2

% Initialize log acceptance probability to -Inf (automatic rejection)
log_alpMH = -inf(n, 1);

% Identify stationary candidates
is_stationary = abs(phic) < STATIONARITY_BOUND;

if any(is_stationary)
    idx = is_stationary;
    
    % Calculate g(phic) - g(phi) only for stationary candidates
    
    % Numerical Stability: Use log1p(-x) instead of log(1-x) for better precision near x=1
    
    % Term 1: 0.5 * [log(1-phic^2) - log(1-phi^2)]
    log_det_diff = 0.5 * (log1p(-phic(idx).^2) - log1p(-phi(idx).^2));
    
    % Term 2: -0.5 * h1^2/sig2 * [(1-phic^2) - (1-phi^2)]
    % Simplified: -0.5 * h1^2/sig2 * [phi^2 - phic^2]
    quad_diff = -0.5 * (h1_sq(idx) .* iSig2(idx)) .* (phi(idx).^2 - phic(idx).^2);
    
    log_alpMH(idx) = log_det_diff + quad_diff;
end

% 4. Accept/Reject Proposals (Vectorized)

% Numerical Stability: Perform comparison in log scale
% Accept if log(rand) < log_alpMH
accept = log(rand(n, 1)) < log_alpMH;

% Update phi where accepted
phi(accept) = phic(accept);

% Return flag indicating acceptance (1) or rejection (0)
flag_phi = double(accept);
end