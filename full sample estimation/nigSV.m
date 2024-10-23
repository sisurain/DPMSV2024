% Class for the Normal-Inverse-Gamma (NIG) distribution, utilized in the 
% Dirichlet Process with Stochastic Volatility (SV).
%
% In addition to the primary data inputs (X), the SV component (h) and 
% the cross-sectional covariance matrix (Sig) from the VAR model are
% needed.

classdef nigSV
        properties
         
         nu_
         iV_ % we track inverse of V
         mu_
         U_ % we track U, instead of S
     end
     
     methods
         function obj = nigSV(mu,V,nu,S) % hyperparameters
             iV = V\speye(size(V,1));
             U = S+mu'*iV*mu/2;
             obj.iV_ = iV;
             obj.mu_ = mu;
             obj.nu_ = nu;
             obj.U_ = U;
         end
         
         function obj = clone(obj)
         end
         
         function obj = addData(obj, X, h, Sig) % X is n by d, h is n by 1, Sig is d by d
             iV0 = obj.iV_;
             mu0 = obj.mu_;
             nu0 = obj.nu_;
             U0 = obj.U_;
             
             exph = exp(-h/2);
             exph2 = exp(-h);
             n = size(X,1);
             d = size(X,2);
             
             nu = nu0+n*d/2;
             iV = iV0+sum(exph2)*(Sig\speye(d));
             mu = iV\(Sig\(sum(exph2.*X,1)') + iV0*mu0);
             %{
             R = chol(Sig)*((exph.*X)');
             r = sum(dot(R,R,1));
             U = U0+r/2;
             %}
             U = U0;
             for i = 1:size(h)
                 U = U+(exph(i)*X(i,:))*(Sig\((exph(i)*X(i,:))'))/2;
             end
             
             obj.iV_ = iV;
             obj.mu_ = mu;
             obj.nu_ = nu;
             obj.U_ = U;
         end
        
         function obj = addSample(obj, x, h, Sig) % x is 1 by d
             iV = obj.iV_;
             mu = obj.mu_;
             nu = obj.nu_;
             U = obj.U_;
             
             exph = exp(-h/2);
             exph2 = exp(-h);
             d = size(x,2);
             
             iV = iV+exph2*(Sig\speye(d));
             mu = iV\((iV-exph2*(Sig\speye(d)))*mu+Sig\(exph2*x'));
             nu = nu+d/2;
             U = U+(exph*x)*(Sig\((exph*x)'))/2;
             
             obj.iV_ = iV;
             obj.mu_ = mu;
             obj.nu_ = nu;
             obj.U_ = U;
         end
         
         function obj = delSample(obj, x, h, Sig)
             iV = obj.iV_;
             mu = obj.mu_;
             nu = obj.nu_;
             U = obj.U_;
             
             exph = exp(-h/2);
             exph2 = exp(-h);
             d = size(x,2);
             
             iV = iV-exph2*(Sig\speye(d));
             mu = iV\((iV+exph2*(Sig\speye(d)))*mu-Sig\(exph2*x'));
             nu = nu-d/2;
             U = U-(exph*x)*(Sig\((exph*x)'))/2;
             
             obj.iV_ = iV;
             obj.mu_ = mu;
             obj.nu_ = nu;
             obj.U_ = U;
         end
         
         function y = logPredPdf(obj, X, h, Sig) 
             % predictive densiity is multivariate T(exph.*mu, ((U-mu'*iV*mu/2)/nu)*(eye(d)+exph2*Sig\V)*Sig)
             % v is degrees of freedom = 2*nu 
             iV = obj.iV_;
             mu = obj.mu_;
             nu = obj.nu_;
             U = obj.U_;
             
             exph = exp(-h/2);
             exph2 = exp(-h);
             d = size(X,2);
             v = 2*nu;
             U = sqrt((U-mu'*iV*mu/2)/nu)*chol((eye(d)+exph2*((iV*Sig)\eye(d)))*Sig);
             
             X = bsxfun(@minus,exph.*X',exph.*mu);
             Q = U'\X;
             q = dot(Q,Q,1); % quadratic term (M distance)
             o = -log(1+q/v)*((v+d)/2);
             c = gammaln((v+d)/2)-gammaln(v/2)-(d*log(v*pi)+2*sum(log(diag(U))))/2;
             y = c+o;
         end
     end
end
