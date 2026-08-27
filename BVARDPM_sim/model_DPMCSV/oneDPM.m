function [z, Theta, nk] = oneDPM(X, h, Sig, alpha, theta)
% One time (for initialization) collapsed Gibbs sampling for Dirichlet Process (infinite) mixture model. 
%
% Input: 
%   X: T x d data matrix
%   h: T x 1 vector of log-volatilities
%   Sig : d x d cross-sectional covariance matrix
%   alpha: concentration parameter for Dirichlet Process prior
%   theta: class object for prior (e.g., nigSV handle class)
% Output:
%   z: 1 x T cluster label
%   Theta: 1 x K cell array of trained components (nigSV objects)
%   nk: 1 x K component count vector

T = size(X,1);
Theta = {}; % Initialize empty cell array for cluster parameters
nk = [];    % Initialize empty vector for cluster counts
z = zeros(1,T); % Initialize assignments

% Sequential scan in random order
for i = randperm(T)
    x = X(i,:);
    ht = h(i);
    
    % Calculate log probabilities for existing clusters
    
    % Robustness check: Ensure Theta/nk are not empty before calculating Pk
    if ~isempty(nk)
        Pk = log(nk) + cellfun(@(t) t.logPredPdf(x,ht,Sig), Theta);
    else
        % Handle the first iteration robustly
        Pk = [];
    end
        
    % Calculate log probability for a new cluster
    P0 = log(alpha) + theta.logPredPdf(x,ht,Sig);
    
    p = [Pk, P0];
    
    % Normalize probabilities (Log-Sum-Exp trick for stability)
    p = p - max(p,[],2);
    p = exp(p);
    p = p./sum(p,2);
    
    % Sample new assignment
    kk = randmn(p);
    
    % Update cluster parameters and counts
    if kk == numel(Theta)+1
        % New cluster: Clone the prior (essential for Handle class)
        Theta{kk} = theta.clone();
        % Add the sample (modifies the clone in place)
        Theta{kk}.addSample(x,ht,Sig); 
        nk = [nk,1];
    else
        % Existing cluster: Add the sample in place (efficient for Handle class)
        Theta{kk}.addSample(x,ht,Sig);
        nk(kk) = nk(kk)+1;
    end
    z(i) = kk;
end
end