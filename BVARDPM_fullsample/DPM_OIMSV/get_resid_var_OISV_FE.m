% get_resid_var_OISV_FE.m (Improved Version)
function sig2 = get_resid_var_OISV_FE(Y0, Y, p_ar)
% Computes residual variances of univariate AR(p) models WITHOUT intercept.
% Features robust estimation using Ridge fallback and data-driven variance fallback.
% Inputs:
%   Y0: Pre-sample data (T0 x n)
%   Y: Estimation sample data (T x n)
%   p_ar: AR lag length (optional, default=4)

% 1. Configuration and Validation
if nargin < 3
    p_ar = 4; % Default AR lag length
end

[T, n] = size(Y);
T0 = size(Y0, 1);
sig2 = zeros(n, 1);

if T0 < p_ar
    error('get_resid_var_OISV_FE:NotEnoughData', ...
          'Not enough initial data Y0 (%d rows) for AR(%d) estimation.', T0, p_ar);
end

% Ensure positive degrees of freedom
if T <= p_ar
     error('get_resid_var_OISV_FE:NotEnoughData', ...
          'Estimation sample T (%d obs) must be larger than AR lag length (%d).', T, p_ar);
end

% 2. Data Preparation
% Combine data: Use only the necessary part of Y0 (last p_ar rows)
tmpY = [Y0(end-p_ar+1:end, :); Y];
T_total = size(tmpY, 1); % T + p_ar

% Define stability thresholds
MIN_VAR = 1e-12;
RCOND_THRESHOLD = 1e-15;

% 3. Estimation Loop
for i = 1:n
    % Define the dependent variable Yi (T x 1)
    Yi = tmpY(p_ar+1:end, i);
    
    % Construct the regressor matrix Xi (T x p_ar)
    Xi = zeros(T, p_ar);
    for j = 1:p_ar
        % Lag j: Indices from (p_ar - j + 1) to (T_total - j)
        Xi(:, j) = tmpY(p_ar-j+1 : T_total-j, i);
    end
    
    % OLS estimation without intercept
    % Robustness: Check conditioning and use Ridge fallback if necessary
    XtX = Xi'*Xi;
    if rcond(XtX) < RCOND_THRESHOLD
        % Use slight ridge regularization if highly collinear
        lambda_ridge = 1e-8;
        beta = (XtX + lambda_ridge*eye(p_ar)) \ (Xi'*Yi);
    else
        % Use standard backslash operator
        beta = Xi \ Yi;
    end
    
    resid = Yi - Xi*beta;
    
    % Calculate variance (T-p_ar denominator for unbiased estimate)
    variance_i = sum(resid.^2) / (T - p_ar);
    
    % Robustness: Enhanced fallback strategy
    if variance_i <= MIN_VAR || isnan(variance_i) || ~isreal(variance_i)
        % Fallback 1: Use the unconditional variance of the series in the estimation sample
        uncond_var = var(Y(:,i));
        if uncond_var > MIN_VAR && ~isnan(uncond_var) && isreal(uncond_var)
            sig2(i) = uncond_var;
        else
            % Fallback 2: Absolute constant if series is constant or data is invalid
            sig2(i) = 1e-6;
        end
    else
        sig2(i) = variance_i;
    end
end
end