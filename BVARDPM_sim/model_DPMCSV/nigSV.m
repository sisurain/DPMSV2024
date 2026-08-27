% Class for the Normal-Inverse-Gamma (NIG) distribution, utilized in the 
% Dirichlet Process with Stochastic Volatility (SV).
%
% Inherits from 'handle' for efficiency (modification in place) and 
% 'matlab.mixin.Copyable' for proper cloning.
%
% Parameters tracked:
% nu (Shape), iV (Precision V^-1), mu (Mean), U (Transformed Scale S + 0.5*mu'*iV*mu)

classdef nigSV < handle & matlab.mixin.Copyable
    properties
        nu_  % Shape parameter (nu)
        iV_  % Precision matrix (inverse of V)
        mu_  % Mean parameter (mu)
        U_   % Transformed scale parameter (U)
    end
    
    methods
        %% Constructor
        function obj = nigSV(mu, V, nu, S)
            % Initialize the object with hyperparameters
            d = size(V, 1);
            % Calculate precision efficiently using backslash and eye() for dense matrix
            iV = V \ eye(d); 
            % Calculate transformed scale U
            U = S + mu' * iV * mu / 2;
            
            obj.iV_ = iV;
            obj.mu_ = mu;
            obj.nu_ = nu;
            obj.U_ = U;
        end
        
        %% Clone (Copy) Method
        function newObj = clone(obj)
            % Creates a deep copy of the object.
            % Relies on the protected 'copy' method from matlab.mixin.Copyable
            newObj = copy(obj);
        end
        
        %% Batch Update (Vectorized)
        function obj = addData(obj, X, h, Sig)
            % Updates the NIG parameters by adding a batch of data.
            % X: N_batch x d matrix of data
            % h: N_batch x 1 vector of log-volatilities
            
            [N_batch, d] = size(X);
            iSig = Sig \ eye(d);
            
            % Calculate weights
            exph = exp(-h/2); % b_i
            exph2 = exp(-h);  % c_i
            
            % 1. Update U (Vectorized quadratic form calculation)
            X_scaled = exph .* X; % b_i * X_i
            
            % Efficiently calculate the sum of quadratic forms.
            % We use the Cholesky approach for stability and speed.
            try
                % X_std = X_scaled * inv(chol(Sig, 'lower')')
                % Use the forward slash operator / for efficient solving
                X_std = X_scaled / chol(Sig, 'lower')'; 
                % Sum of squares gives the sum of quadratic forms
                % sum(..., 'all') requires recent MATLAB. For compatibility:
                quad_form_sum = sum(sum(X_std.^2));
            catch
                % Fallback if Sig is numerically unstable (use iSig calculated earlier)
                quad_form_sum = sum(sum((X_scaled * iSig) .* X_scaled));
            end
            obj.U_ = obj.U_ + 0.5 * quad_form_sum;

            % 2. Calculate terms for mu update (uses OLD iV and OLD mu)
            % Efficient calculation of Sum(c_i * X_i') using matrix mult
            weighted_data_sum_X = (exph2' * X)';
            % iSig * Sum(c_i * X_i') + iV0*mu0
            mu_update_term = iSig * weighted_data_sum_X + obj.iV_ * obj.mu_;

            % 3. Update iV
            % iV = iV0 + Sum(c_i) * iSig
            obj.iV_ = obj.iV_ + sum(exph2) * iSig;
            
            % 4. Update mu (using NEW iV)
            obj.mu_ = obj.iV_ \ mu_update_term;

            % 5. Update nu
            obj.nu_ = obj.nu_ + N_batch * d / 2;
        end
        
        %% Sequential Update (Add single sample)
        function obj = addSample(obj, x, h, Sig)
            % Updates the NIG parameters by adding a single sample.
            % Optimized for Handle class (modification in place).
            
            d = size(x, 2);
            iSig = Sig \ eye(d);
            
            % Calculate weights
            exph = exp(-h/2); % b
            exph2 = exp(-h);  % c
            W = exph2 * iSig; % Weighted precision contribution c*iSig
            
            % 1. Update U
            X_scaled = exph * x;
            % Quadratic form: (b*x) * iSig * (b*x)'
            quad_form = X_scaled * iSig * X_scaled';
            obj.U_ = obj.U_ + quad_form / 2;
            
            % 2. Calculate terms for mu update (uses OLD iV and OLD mu)
            % W*x' + iV0*mu0
            mu_update_term = W*x' + obj.iV_ * obj.mu_;

            % 3. Update iV
            obj.iV_ = obj.iV_ + W;
            
            % 4. Update mu (using NEW iV)
            obj.mu_ = obj.iV_ \ mu_update_term;

            % 5. Update nu
            obj.nu_ = obj.nu_ + d / 2;
        end
        
        %% Sequential Downdate (Remove single sample)
        function obj = delSample(obj, x, h, Sig)
            % Updates the NIG parameters by removing a single sample.

            d = size(x, 2);
            iSig = Sig \ eye(d);
            
            % Calculate weights
            exph = exp(-h/2);
            exph2 = exp(-h);
            W = exph2 * iSig;
            
            % 1. Update U
            X_scaled = exph * x;
            quad_form = X_scaled * iSig * X_scaled';
            obj.U_ = obj.U_ - quad_form / 2;
            
            % 2. Calculate terms for mu update (uses OLD iV and OLD mu)
            % iV0*mu0 - W*x'
            mu_update_term = obj.iV_ * obj.mu_ - W*x';

            % 3. Update iV
            obj.iV_ = obj.iV_ - W;
            
            % 4. Update mu (using NEW iV)
            obj.mu_ = obj.iV_ \ mu_update_term;

            % 5. Update nu
            obj.nu_ = obj.nu_ - d / 2;
        end
        
        %% Log Posterior Predictive Density (Multivariate Student's t)
        function y = logPredPdf(obj, X, h, Sig)
            % Calculates the log predictive density for the scaled observation b*X.
            % Distribution: t_v(b*mu, ScaleMatrix), where v = 2*nu.
            % ScaleMatrix = (S/nu) * (Sig + c_new*V)
            
            % Retrieve properties
            iV = obj.iV_; mu = obj.mu_; nu = obj.nu_; U_prop = obj.U_;
            
            % Calculate weights
            exph = exp(-h/2); % b_new
            exph2 = exp(-h);  % c_new
            
            d = size(X, 2);
            v = 2*nu; % Degrees of freedom

            % 1. Recover standard scale S from transformed scale U
            % S = U - 0.5 * mu'*iV*mu
            S = U_prop - 0.5 * (mu' * iV * mu);
            
            % Robustness check: S must be positive
            if S <= 0
                y = -inf; return;
            end
             
            % 2. Calculate V = inv(iV)
            V = iV \ eye(d);

            % 3. Calculate the Correct Scale Matrix
            ScaleMatrix = (S/nu) * (Sig + exph2 * V);
             
            % 4. Calculate Cholesky factor L (Robustly)
            [L, flag] = chol(ScaleMatrix, 'lower');
            if flag ~= 0
                % Handle cases where numerical instability makes the matrix non-PD
                y = -inf; 
                return; 
            end

            % 5. Calculate Mahalanobis distance
            % Difference between scaled data and scaled mean
            % X_diff = (b*X)' - b*mu
            % Assumes X is a single observation (1xd) as used in the DPM loop
            X_diff = (exph .* X)' - exph .* mu;
            
            % Solve L*Q = X_diff for Q
            Q = L \ X_diff; 
            
            % q = Q'*Q (Mahalanobis distance squared)
            q = Q' * Q; 
             
            % 6. Calculate log density
            % Log quadratic term
            o = -log(1 + q/v) * ((v+d)/2);
            
            % Log normalizing constant
            % Calculate log determinant of ScaleMatrix using L: 2*sum(log(diag(L)))
            log_det = 2 * sum(log(diag(L)));
            c = gammaln((v+d)/2) - gammaln(v/2) - (d*log(v*pi) + log_det)/2;
            
            y = c + o;
        end
    end
end