% This script samples the common stochastic volatility h using an 
% Acceptance-Rejection Metropolis-Hastings (AR-MH) algorithm with 
% a Gaussian proposal constructed via Newton-Raphson (IWLS).
%
% References: Chan (2020, JBES).

function [h,is_accept] = sample_h(s2,rho,sigh2,h,n)
T = size(s2,1);
is_accept = 0;
max_iter = 50; % Maximum iterations for Newton-Raphson
tol = 1e-5;    % Convergence tolerance (Increased precision)

% 1. Construct Prior Precision Matrix (HiSH)
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
% Diagonal elements of the precision of the innovations (Efficiently using spdiags)
prior_diag = [(1-rho^2)/sigh2; 1/sigh2*ones(T-1,1)];
% Use spdiags for efficient construction of the diagonal matrix
Dinv = spdiags(prior_diag, 0, T, T);
HiSH = Hrho'*Dinv*Hrho;

% 2. Newton-Raphson Iteration (IWLS)
errh = inf; 
ht = h; % Initialize at the previous state
iter = 0;

while errh > tol && iter < max_iter
    iter = iter + 1;
    eht = exp(ht);
    % Add small epsilon for numerical stability if eht is near zero
    sieht = s2./(eht + 1e-8); 
    
    % Gradient (fh) and Negative Hessian (Gh) of the log-likelihood
    fh = -n/2 + .5*sieht;
    Gh = .5*sieht;
    
    % Posterior Precision (Efficiently using spdiags)
    Kh = HiSH + spdiags(Gh, 0, T, T);
    
    % Newton step (Solving the linear system)
    newht = Kh\(fh+Gh.*ht);
    
    errh = max(abs(newht-ht));
    ht = newht;          
end

% Check for stability issues before Cholesky
[CKh, flag_chol] = chol(Kh,'lower');
if flag_chol ~= 0
    % If Kh is not positive definite, the approximation failed.
    % Return the previous h to maintain the Markov chain validity.
    return; 
end

% 3. AR-MH Step

% Define the log-target density f(h) (ignoring constants)
% f(h) = log p(h) + log p(s2|h) 
%      = -0.5*h'*HiSH*h - n/2*sum(h) - 0.5*exp(-h)'*s2
log_target = @(x) -.5*x'*HiSH*x - n/2*sum(x) - .5*(exp(-x)'*s2);

% Set the AR bound logc (using the mode found by NR)
% Use an inflation factor (e.g., log(3)) to increase AR acceptance rate
inflation_factor = log(3); 
logc = log_target(ht) + inflation_factor; 

% AR-step: Draw from the proposal until accepted
flag_AR = 0;
while flag_AR == 0
    % Draw candidate hc from proposal N(ht, Kh^-1)
    hc = ht + CKh'\randn(T,1);
    
    % Calculate R(h) = log_target(h) - log_proposal(h)
    % log_proposal(h) = -0.5*(h-ht)'*Kh*(h-ht)
    
    % Efficient calculation of the quadratic term for the proposal
    diff_hc = hc - ht;
    % We can reuse the Cholesky factor: (diff'*Kh*diff) = (CKh*diff)'*(CKh*diff)
    tmp_hc = CKh * diff_hc;
    log_R_hc = log_target(hc) - (-.5*(tmp_hc'*tmp_hc));
    
    % AR acceptance ratio
    alpARc = log_R_hc - logc;
    
    if alpARc > log(rand)
        flag_AR = 1;
    end
end        

% MH-step: Ensure detailed balance
% Calculate R(h) for the current state h
diff_h = h - ht;
tmp_h = CKh * diff_h;
log_R_h = log_target(h) - (-.5*(tmp_h'*tmp_h));

% MH acceptance ratio (log scale): R(hc) - R(h)
% CORRECTION: This replaces the flawed logic in the original script.
alpMH = log_R_hc - log_R_h;

% Accept or reject
if alpMH > log(rand)
    h = hc; 
    is_accept = 1;
end
end