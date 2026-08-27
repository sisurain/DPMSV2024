% construct_Sigt_DPM_OISV.m (Improved Version)
function Sigt = construct_Sigt_DPM_OISV(H, B0, lam)
% Constructs the full time-varying covariance matrix lambda_t * Sigma_t
% Inputs:
%   H: T x n matrix of log-volatilities
%   B0: n x n contemporaneous matrix
%   lam: T x 1 vector of DPM scaling factors (lambda_t)

% 1. Input Validation and Initialization
[T, n] = size(H);

assert(size(B0, 1) == n && size(B0, 2) == n, 'construct_Sigt_DPM_OISV:DimensionMismatch', ...
    'B0 must be a square matrix (n x n) with dimensions matching the columns of H.');
assert(isvector(lam) && length(lam) == T, 'construct_Sigt_DPM_OISV:DimensionMismatch', ...
    'lam must be a vector with length matching the rows of H (T).');

% Ensure lam is a column vector for consistency
if isrow(lam)
    lam = lam';
end

Sigt = zeros(n, n, T);

% 2. Robust Inversion of B0
% Use the backslash operator (mldivide) for efficient inversion: iB0*B0 = I
try
    iB0 = B0 \ eye(n);
catch ME
     error('construct_Sigt_DPM_OISV:InversionFailed', ...
          'B0 matrix is singular and cannot be inverted. Error: %s', ME.message);
end

% Check reciprocal condition number for stability diagnostics
r_cond = rcond(B0);
if r_cond < 1e-15
    % Warn if B0 is nearly singular, as iB0 will be inaccurate.
    warning('construct_Sigt_DPM_OISV:SingularMatrix', ...
            'B0 matrix is close to singular (RCOND = %e). Covariance results may be inaccurate.', r_cond);
end

% 3. Numerical Stability and Optimization
H_MAX = 50; % Define bound to prevent overflow (exp(50) is very large)
H_clamped = H;
H_clamped(H_clamped > H_MAX) = H_MAX;
% Note: Underflow (H very small negative) is generally fine as exp(H) -> 0.

% Pre-calculate exponentiated volatilities
expH = exp(H_clamped);

% 4. Calculation Loop
for t=1:T
    % Optimization: Combine DPM scaling (lambda_t) and SV scaling (D_t)
    % Dt_full_diag = lam(t) * exp(H(t,:)) (The diagonal of lambda_t * Dt)
    Dt_full_diag = lam(t) * expH(t,:);
    
    % Calculate Full Covariance using broadcasting: (iB0 .* Dt_full_diag) * iB0'
    % This efficiently calculates B0^{-1} * (lambda_t * Dt) * (B0^{-1})'
    S_t = (iB0 .* Dt_full_diag) * iB0';
    
    % Robustness: Enforce Symmetry (addresses floating-point inaccuracies)
    % Ensures the matrix is perfectly symmetric
    Sigt(:,:,t) = 0.5 * (S_t + S_t');
end
end