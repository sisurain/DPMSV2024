% File: oneDPM_MSV.m
function [z, Theta, nk] = oneDPM_MSV(X, B, H, alpha, theta)
% One time (for initialization) collapsed Gibbs sampling for DPM-MSV.
% Input: 
%   X: T x n data matrix
%   B: n x n contemporaneous matrix
%   H: T x n matrix of log-volatilities
%   alpha: concentration parameter
%   theta: class object for prior (nigMSV handle class)

T = size(X,1);
Theta = {};
nk = [];
z = zeros(1,T);

% Sequential scan in random order
for i = randperm(T)
    x = X(i,:);
    h_vec = H(i,:); % Extract the volatility vector for time i
    
    % Calculate log probabilities for existing clusters
    if ~isempty(nk)
        % Pass B and the specific volatility vector h_vec
        Pk = log(nk) + cellfun(@(t) t.logPredPdf(x, B, h_vec), Theta);
    else
        Pk = [];
    end
        
    % Calculate log probability for a new cluster
    P0 = log(alpha) + theta.logPredPdf(x, B, h_vec);
    
    p = [Pk, P0];
    
    % Normalize and Sample
    p = p - max(p,[],2);
    p = exp(p);
    p = p./sum(p,2);
    kk = randmn(p); % Assumes randmn.m is available
    
    % Update cluster parameters and counts
    if kk == numel(Theta)+1
        % New cluster: Clone the prior
        Theta{kk} = theta.clone();
        Theta{kk}.addSample(x, B, h_vec); 
        nk = [nk,1];
    else
        % Existing cluster: Add the sample in place
        Theta{kk}.addSample(x, B, h_vec);
        nk(kk) = nk(kk)+1;
    end
    z(i) = kk;
end
end