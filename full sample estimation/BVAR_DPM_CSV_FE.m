% This script estimates the BVAR-DPM-CSV model 
% with FRED datasets ending in 2022-12, excluding the 'a0' parameter.

clear; clc; close all;

p = 4; % VAR lag length. If p > 4, adjust Y0 and shortY accordingly.
%p = 12; 
nsims = 1000; % Number of MCMC simulations.
burnin = 500; % Burn-in period.

%% Load Data
data = readmatrix('2022-12-ccmm16-r.txt'); % Read data from file.

Y0 = data(1:4,:); % Initial conditions using the first 4 observations.
shortY = data(5:end,:); % Dataset after initial observations.
%Y0 = data(1:12,:); % Alternative initial condition with 12 lags.
%shortY = data(13:end,:);
[T,n] = size(shortY); % T: Time periods, n: Variables.   
k = n*p;      

%% Set Prior Parameters
S0 = eye(n); nu0 = n+3;
nuh0 = 1e4; Sh0 = .001*(nuh0-1); 
%nuh0 = 0.01; Sh0 = .1*(nuh0-1);
rho0 = .9; Vrho = .2^2; % Rho prior.
nuub = 100; % Upperbound for nu.
construct_prior_A_woa0; % Build prior without 'a0'.

a0 = 20; b0 = 8; % Hyperparameters for alpha (DPM).

mu0 = zeros(n,1); V0 = 100*eye(n); 
%nu0_dp = 1; S0_dp = 10;
nu0_dp = 10; S0_dp = 100;
prior_dp = nigSV(mu0,V0,nu0_dp,S0_dp);

%% Construct Matrix X
X = zeros(T,n*p); 
for i=1:p
    X(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end

%% Initialize Storage Matrices for MCMC
store_Sig = zeros(n,n); 
store_A = zeros(k,n);
store_h = zeros(nsims,T);
store_lam = zeros(nsims,T);
store_mu = zeros(T,n,nsims);
store_z = zeros(nsims,T);
store_theta = zeros(nsims,4);

store_res1 = zeros(T,n,nsims);
store_res2 = zeros(T,n,nsims);

store_yt1 = zeros(nsims,n);

counth = 0; countrho = 0; countnu = 0;
nugrid = linspace(2,nuub,700)';
store_pnu = zeros(700,1);

%% Initialize the Chain
h = zeros(T,1);
nu = 5;
rho = .8;
sigh2 = .01;
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
lam = 1./gamrnd(nu/2,2/nu,T,1);
mu = zeros(T,n);

alpha = 1e5; 

[z,Theta,nk] = oneDPM(shortY,h,S0,alpha,prior_dp); % Initial cluster assignments.

%% MCMC Starts Here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000)); % Set random seed.
disp('Starting MCMC for BVAR-DPM-CSV....');
start_time = clock;

for isim = 1:nsims + burnin
  
    %% sample Sig and A    
    iOm = sparse(1:T,1:T,exp(-h)./lam);
    XiOm = X'*iOm;
    KA = sparse(1:k,1:k,1./VA0) + XiOm*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOm*(shortY-mu));
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + (shortY-mu)'*iOm*(shortY-mu) ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);    
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig'; 
    
    %% sample h
    U = shortY - X*A -mu;
    tmp = (U/CSig');
    s2 = sum(tmp.^2,2)./lam;
    [h,flag] = sample_h(s2,rho,sigh2,h,n);
    counth = counth + flag;
    
    %% sample lam and mu
    X_dp = shortY - X*A;
    for j = 1:length(Theta)
        idxT = find(z==j);
        Theta{j} = prior_dp.addData(X_dp(idxT,:),h(idxT),Sig);
    end
   
    for i = randperm(T)
        x = X_dp(i,:);
        hi = h(i);
        kk = z(i);
        Theta{kk} = Theta{kk}.delSample(x,hi,Sig);
        nk(kk) = nk(kk)-1;

        Pk = log(nk)+cellfun(@(t) t.logPredPdf(x,hi,Sig), Theta);
        P0 = log(alpha)+prior_dp.logPredPdf(x,hi,Sig);
        pp = [Pk,P0];
        pp = pp - max(pp,[],2);
        pp = exp(pp);
        pp = pp./sum(pp,2);
    
        kk = randmn(pp);
        if kk == numel(Theta)+1 % add extra cluster
            Theta{kk} = prior_dp.clone.addSample(x,hi,Sig);
            nk = [nk,1];
        else
            Theta{kk} = Theta{kk}.addSample(x,hi,Sig);
            nk(kk) = nk(kk)+1;
        end
        z(i) = kk; % update z
        
        idx0 = find(nk==0); % remove empty cluster

        if ~isempty(idx0)
            Theta(idx0) = [];
            nk(idx0) = [];
            which = z>idx0;
            z(which) = z(which)-1;
        end     
    end
    
    K = length(Theta);
    
    lamK = zeros(K,1); muK = zeros(K,n);
    for i = 1:K
        Ui = Theta{i}.U_;
        iVi = Theta{i}.iV_;
        mui = Theta{i}.mu_;
        nui = Theta{i}.nu_;
        S = Ui - mui'*iVi*mui/2;
        lamK(i) = 1./gamrnd(nui/2,2/S,1,1);
        muK(i,:) = (mui + chol(iVi/lamK(i),'lower')\randn(n,1))';
    end
    
    zK = sparse([1:T]',z',ones(T,1),T,K);   
    lam = zK*lamK; mu = zK*muK;

    %% sample sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];    
    sigh2 = 1/gamrnd(nuh0+T/2,1/(Sh0 + sum(eh.^2)/2));

    %% sample rho
    Krho = 1/Vrho + sum(h(1:T-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:T-1)'*h(2:T)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.9999
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;
            countrho = countrho+1;
            Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
        end
    end 
     
    %% sample alpha

    xi = beta(alpha+1,T);
    pi = (a0+K-1)/((a0+K-1)+T*(b0-log(xi)));
    if rand < pi
        alpha = gamrnd(a0+K, 1/(b0-log(xi)));
    else
        alpha = gamrnd(a0+K-1, 1/(b0-log(xi)));
    end
    
    if isim > burnin
        isave = isim - burnin; 
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;
        store_h(isave,:) = h';
        store_lam(isave,:) = lam';
        store_theta(isave,:) = [alpha K rho sigh2];
        store_z(isave,:) = z;
        store_mu(:,:,isave) = mu;      
    end
    
    if ( mod(isim, 100) ==0 )
        disp(  [ num2str(isim) ' loops... ' ] )
    end 
    
end

disp( ['MCMC takes' num2str(etime(clock, start_time)) 'seconds']);
disp(' ');

A_mean = store_A/nsims;
Sig_mean = store_Sig/nsims;
theta_mean = mean(store_theta)';
h_mean = mean(store_h)';
%pnu_mean = store_pnu/nsims;
CSV_mean = mean(exp(store_h/2))';
lam_mean = mean(store_lam)';
%CSVlam = CSV_mean.*lam_mean;

figure;
colormap('hsv');
imagesc(A_mean);
colorbar;
box off;
title('Heat map of the VAR coefficients');    

T_id = linspace(1959+7/12,2022+11/12,T)';

figure; 
plot(T_id, CSV_mean); box off;
title('Posterior mean of $e^{h/2}$','Interpreter','latex');
box off; xlim([T_id(1)-1 T_id(end)+1]);
xticks(1955:5:2020)

figure; 
plot(T_id, lam_mean); box off;
title('Posterior mean of $\lambda$','Interpreter','latex');
box off; xlim([T_id(1)-1 T_id(end)+1]);
xticks(1955:5:2020)

lam_mean1 = lam_mean;
lam_mean1(lam_mean1>20) = 20;
find(lam_mean1==20)

figure; 
plot(T_id, lam_mean1); box off;
title('Posterior mean of $\lambda$ with values less than 20','Interpreter','latex');
box off; xlim([T_id(1)-1 T_id(end)+1]);
xticks(1955:5:2020)

[unique_z,~,ix] = unique(store_z(end,:));
C1 = accumarray(ix,1).'
for i = 2:length(unique_z)
    find(store_z(end,:)==i)
end

ans1 = [min(lam_mean), max(lam_mean), median(lam_mean)]
ans2 = find(lam_mean>20)
ans3 = find(lam_mean>100)
ans4 = sum(lam_mean<2*median(lam_mean))

