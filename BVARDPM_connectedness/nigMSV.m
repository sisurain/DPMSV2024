% File: nigMSV.m
% Class for the Normal-Inverse-Gamma (NIG) distribution, adapted for DPM-MSV.
% Inputs use the Cholesky structure (B, H) where Sigma_t^-1 = B' * diag(exp(-H_t)) * B.

classdef nigMSV < handle & matlab.mixin.Copyable
    properties
        nu_
        iV_
        mu_
        U_
    end
    
    methods
        %% Constructor
        function obj = nigMSV(mu, V, nu, S)
            d = size(V, 1);
            iV = V \ eye(d); 
            U = S + mu' * iV * mu / 2;
            obj.iV_ = iV;
            obj.mu_ = mu;
            obj.nu_ = nu;
            obj.U_ = U;
        end
        
        %% Clone (Copy) Method
        function newObj = clone(obj)
            newObj = copy(obj);
        end
        
        %% Batch Update (Efficiently using B and H)
        function obj = addData(obj, X, B, H)
            % X: N_batch x n data matrix
            % B: n x n contemporaneous matrix
            % H: N_batch x n matrix of log-volatilities

            [N_batch, n] = size(X);
            iH_weights = exp(-H); 

            % 1. Update U (Vectorized Quadratic Form)
            % Orthogonalize residuals: U_orth = (X * B')' (n x N_batch)
            U_orth = (X * B')'; 
            % sum( sum( weights .* (U_orth.^2) ) )
            quad_form_sum = sum(sum( iH_weights' .* (U_orth.^2) ));
            obj.U_ = obj.U_ + 0.5 * quad_form_sum;

            % 2. Calculate Sum_t(Sigma_t^-1) and Sum_t(Sigma_t^-1 * X_t')
            % Optimized loop leveraging the structure iSig_t = (B' * D * B).
            Sum_iSig = zeros(n, n);
            Sum_iSigX = zeros(n, 1);

            for t = 1:N_batch
                % Efficient calculation of iSig_t = (B' .* w_t) * B
                iSig_t = (B' .* iH_weights(t,:)) * B;
                
                Sum_iSig = Sum_iSig + iSig_t;
                Sum_iSigX = Sum_iSigX + iSig_t * X(t,:)';
            end

            % 3. Calculate mu update term (uses OLD iV and OLD mu)
            mu_update_term = Sum_iSigX + obj.iV_ * obj.mu_;

            % 4. Update iV
            obj.iV_ = obj.iV_ + Sum_iSig;
            
            % 5. Update mu (using NEW iV)
            obj.mu_ = obj.iV_ \ mu_update_term;

            % 6. Update nu
            obj.nu_ = obj.nu_ + N_batch * n / 2;
        end
        
        %% Sequential Update (Add single sample)
        function obj = addSample(obj, x, B, h_vec)
            n = size(x, 2);
            iH_weights = exp(-h_vec);

            % 1. Update U
            U_orth = (x * B')';
            quad_form = sum( iH_weights' .* (U_orth.^2) );
            obj.U_ = obj.U_ + 0.5 * quad_form;
            
            % 2. Calculate Sigma_t^-1
            iSig_t = (B' .* iH_weights) * B;

            % 3. Calculate mu update term
            mu_update_term = iSig_t * x' + obj.iV_ * obj.mu_;

            % 4. Update iV
            obj.iV_ = obj.iV_ + iSig_t;
            
            % 5. Update mu
            obj.mu_ = obj.iV_ \ mu_update_term;

            % 6. Update nu
            obj.nu_ = obj.nu_ + n / 2;
        end
        
        %% Sequential Downdate (Remove single sample)
        function obj = delSample(obj, x, B, h_vec)
            n = size(x, 2);
            iH_weights = exp(-h_vec);

            % 1. Update U
            U_orth = (x * B')';
            quad_form = sum( iH_weights' .* (U_orth.^2) );
            obj.U_ = obj.U_ - 0.5 * quad_form;
            
            % 2. Calculate Sigma_t^-1
            iSig_t = (B' .* iH_weights) * B;

            % 3. Calculate mu update term
            mu_update_term = obj.iV_ * obj.mu_ - iSig_t * x';

            % 4. Update iV
            obj.iV_ = obj.iV_ - iSig_t;
            
            % 5. Update mu
            obj.mu_ = obj.iV_ \ mu_update_term;

            % 6. Update nu
            obj.nu_ = obj.nu_ - n / 2;
        end
        
        %% Log Posterior Predictive Density (Multivariate Student's t)
        function y = logPredPdf(obj, X, B, H_new)
            % Calculates the log predictive density for X.
            % ScaleMatrix = (S/nu) * (Sigma_new + V)
            
            iV = obj.iV_; mu = obj.mu_; nu = obj.nu_; U_prop = obj.U_;
            n = size(X, 2);
            v = 2*nu;

            % 1. Recover S from U
            S = U_prop - 0.5 * (mu' * iV * mu);
            if S <= 0
                y = -inf; return;
            end
             
            % 2. Calculate V = inv(iV)
            V = iV \ eye(n);

            % 3. Calculate Sigma_new = B^-1 * D_new * B^-1'
            D_new = diag(exp(H_new));
            try
                % Efficient inversion using backslash
                iB = B \ eye(n);
                Sig_new = iB * D_new * iB';
            catch
                % Handle failure if B is singular
                y = -inf; return;
            end

            % 4. Calculate the Correct Scale Matrix
            ScaleMatrix = (S/nu) * (Sig_new + V);
             
            % 5. Calculate Cholesky factor L
            [L, flag] = chol(ScaleMatrix, 'lower');
            if flag ~= 0
                y = -inf; return; 
            end

            % 6. Calculate Mahalanobis distance and log density
            X_diff = X' - mu;
            Q = L \ X_diff; 
            q = Q' * Q; 
             
            o = -log(1 + q/v) * ((v+n)/2);
            log_det = 2 * sum(log(diag(L)));
            c = gammaln((v+n)/2) - gammaln(v/2) - (n*log(v*pi) + log_det)/2;
            
            y = c + o;
        end
    end
end