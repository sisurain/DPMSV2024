% get_C_OISV_FE.m (Vectorized and Robust)
% get_C_OISV_FE.m (Corrected, Vectorized, and Robust)
function [C, idx_kappa1, idx_kappa2] = get_C_OISV_FE(n, p, sig2)
% Constructs the scaling matrix C for the Minnesota-type Horseshoe prior.
% Optimized using vectorization. Adjusted dimensions for k=np (no intercept).

% Input Validation
assert(isvector(sig2) && length(sig2) == n, 'get_C_OISV_FE:InputError', 'sig2 must be a vector of length n.');

k = n*p; % Regressors per equation

% 1. Robustness: Ensure sig2 is strictly positive
MIN_SIG2 = 1e-12;
sig2_robust = sig2;
sig2_robust(sig2_robust < MIN_SIG2) = MIN_SIG2;

% Ensure sig2 is a column vector for consistency
if isrow(sig2_robust)
    sig2_robust = sig2_robust';
end

% 2. Vectorized Lag Scaling (L_scaling)
% Create a vector of lags (k x 1): [1...1, 2...2, ..., p...p] (n repetitions each)
lags = kron((1:p)', ones(n, 1));
L_scaling_vec = 1./(lags.^2); % (k x 1)

% 3. Vectorized Variance Ratios (S)
% Identify the variable index (jj) for each regressor (j)
% vars_jj (k x 1): [1..n, 1..n, ..., 1..n]
vars_jj = kron(ones(p, 1), (1:n)');

% Create matrices of sig2(ii) and sig2(jj) for element-wise division
Sig2_ii = repmat(sig2_robust', k, 1);
Sig2_jj = repmat(sig2_robust(vars_jj), 1, n);

S = Sig2_ii ./ Sig2_jj; % (k x n)

% 4. Identify Own vs Other Lags (Mask_Own)
% Create matrices identifying the equation index (ii) and variable index (jj)
Eq_ii = repmat(1:n, k, 1);
Var_jj = repmat(vars_jj, 1, n);

% Logical mask (k x n): True where equation index equals variable index
Mask_Own = (Eq_ii == Var_jj);

% 5. Construct the Scaling Matrix (C_matrix)
% Initialize using the 'Other lag' rule: L_scaling .* S
% L_scaling_vec (k x 1) broadcasts correctly with S (k x n)
C_matrix = L_scaling_vec .* S;

% Overwrite with the 'Own lag' rule: L_scaling
% CORRECTION: Expand L_scaling into a matrix to match dimensions for indexing
L_scaling_matrix = repmat(L_scaling_vec, 1, n); % (k x n)

C_matrix(Mask_Own) = L_scaling_matrix(Mask_Own);

% 6. Final Output Formatting
C = C_matrix(:); % Vectorize the matrix (Column-major order)
idx_kappa1 = find(Mask_Own(:)); % Find indices for own lags
idx_kappa2 = find(~Mask_Own(:)); % Find indices for other lags
end