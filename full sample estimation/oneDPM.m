function [z, Theta, nk, llh] = oneDPM(X, h, Sig, alpha, theta)
% One time (for initialization) collapsed Gibbs sampling for Dirichlet Process (infinite) mixture model. 
% Input: 
%   X: n x d data matrix
%   h: n x 1 in e^(h/2) of SV
%   Sig : d x d cross-sectional covariance matrix in VAR system
%   alpha: concentration parameter for Dirichlet Process prior
%   theta: class object for prior of component distribution (such as normal inverse Gamma)
% Output:
%   z: 1 x n cluster label
%   Theta: 1 x k structure of trained components (cluster information)
%   w: 1 x k component weight vector
%   llh: loglikelihood
T = size(X,1);
Theta = {};
nk = [];
z = zeros(1,T);
llh = 0;
for i = randperm(T)
    x = X(i,:);
    ht = h(i);
    Pk = log(nk)+cellfun(@(t) t.logPredPdf(x,ht,Sig), Theta);
    P0 = log(alpha)+theta.logPredPdf(x,ht,Sig);
    p = [Pk,P0];
    llh = llh+sum(p-log(T));
    
    p = p - max(p,[],2);
    p = exp(p);
    p = p./sum(p,2);
    
    kk = randmn(p);
    if kk == numel(Theta)+1
        Theta{kk} = theta.clone().addSample(x,ht,Sig);
        nk = [nk,1];
    else
        Theta{kk} = Theta{kk}.addSample(x,ht,Sig);
        nk(kk) = nk(kk)+1;
    end
    z(i) = kk;
end