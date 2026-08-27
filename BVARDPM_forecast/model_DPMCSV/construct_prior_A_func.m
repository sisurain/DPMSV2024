% Constructs the components for a symmetric conjugate Minnesota prior 
% for VAR coefficients A. Can optionally include the intercept 'a0'.
%
% For further reading on the methodology and applications, see:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

function [A0, VA0, sig2] = construct_prior_A_func(Y0, shortY, p, include_intercept)
%
% Inputs:
%   Y0: Pre-sample data (at least p x n)
%   shortY: Estimation sample data (T x n)
%   p: Lag length
%   include_intercept: Boolean (true/false). If true, includes the intercept.
%
% Outputs:
%   A0: Prior mean (k x n)
%   VA0: Prior variances (k x 1). The diagonal of V_prior.
%   sig2: Scaling factors (n x 1)
%
% NOTE: If include_intercept is true, this assumes the intercept '1' is 
% the FIRST column in the main regressor matrix X.

if nargin < 4
    include_intercept = false; % Default to no intercept if not specified
end

[T, n] = size(shortY);

% Define Hyperparameters
c1 = 0.2^2; % Overall tightness for slopes (Standard Minnesota setting)
c2 = 100;   % Prior variance for the intercept (typically large/uninformative)

% Determine the number of regressors (k)
k_slopes = n*p;
if include_intercept
    k = k_slopes + 1;
else
    k = k_slopes;
end

% Initialize
A0 = zeros(k,n); % Prior mean (shrinkage towards zero for all coefficients)
sig2 = zeros(n,1);

% Combine data for AR(p) estimation
if size(Y0, 1) < p
    error('Y0 must contain at least p observations for initialization.');
end
% Use the last p observations from Y0
tmpY = [Y0(end-p+1:end,:); shortY];

%% 1. Calculate scaling factors (sig2) based on AR(p) residuals
% We estimate AR(p) with an intercept to get robust scaling factors, 
% regardless of whether the final VAR includes an intercept.
for i=1:n
    Y_i = tmpY(:,i);
    Y_target = Y_i(p+1:end); % Length T
    
    % Construct lagged regressors dynamically (Generalized for any p)
    X_i = zeros(T, p);
    for lag = 1:p
        X_i(:, lag) = Y_i(p+1-lag : end-lag);
    end
    
    % Add intercept for AR(p) estimation
    Z = [ones(T,1), X_i];
    
    % OLS (efficiently and stably via backslash operator)
    % Solves Z*b = Y_target for b
    tmpb = Z \ Y_target; 
    
    % Calculate residual variance (MSE)
    residuals = Y_target - Z*tmpb;
    sig2(i) = mean(residuals.^2);
end

%% 2. Assign Prior Variances (VA0) for the Slopes (Vectorized)

% Lag vector (lags 1 repeated n times, then 2 repeated n times, etc.)
lags = repelem(1:p, n)';

% Variable index vector (1..n repeated p times)
indices = repmat((1:n)', p, 1);

% Get the corresponding scaling factors
scaling_factors = sig2(indices);

% Calculate VA0 vectorized using harmonic lag decay (l^2)
% Prior Variance = c1 / (l^2 * sigma^2_j)
VA0_slopes = c1 ./ (lags.^2 .* scaling_factors);

%% 3. Assemble the final VA0 vector

if include_intercept
    % If intercept is included, prepend the intercept variance c2
    VA0 = [c2; VA0_slopes];
else
    VA0 = VA0_slopes;
end

end