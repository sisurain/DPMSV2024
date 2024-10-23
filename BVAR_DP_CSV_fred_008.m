% This script estimates the BVAR-DPM-CSV model 
% with FRED datasets ending in 2022-12, excluding the 'a0' parameter.
% Running Forecasts in Parallel


%data = readmatrix(''2022-12-r.txt'');
%data = data(1:203,:);

p = 4; % VAR lag length. If p > 4, adjust Y0 and shortY accordingly.
nsims = 1000;
burnin = 500;

Y0 = data0(1:4,:); % Initial conditions using the first 4 observations.
shortY = data0(5:end,:);
shortY_raw = shortY;

[T,n] = size(shortY);
k = n*p;      

%% Set Prior Parameters
S0 = eye(n); nu0 = n+3;
nuh0 = 1e4; Sh0 = .001*(nuh0-1); 
%nuh0 = 0.01; Sh0 = .1*(nuh0-1);
rho0 = .9; Vrho = .2^2; % Rho prior.
nuub = 100; %% Upperbound for nu
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
%X = [ones(T,1) X]; % with 'a0'.

X_t1 = zeros(1, n*p);
for i=1:p
    X_t1(:,(i-1)*n+1:i*n) = tmpY(end-i+1,:);
end
%X_t1 = [1 X_t1]; % with 'a0'.

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

%% Set Forecast Storage
tH = 24;
nforecast = 100;

store_forecasts = zeros(n,tH,nforecast); 
%tmpyhat1 = zeros(nsims, 2*n);  % [point forecasts, %prelike]
%tmpyhat4 = zeros(nsims, 2*n);
%tmpyhat8 = zeros(nsims, 2*n);
fcstYdraws = zeros(n,tH,nforecast,nsims); 
fcstYhat = zeros(n,tH);
fcstCRPS = zeros(n,tH);

%% MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000)); % Set random seed.
disp('Starting MCMC for BVAR-DPM-CSV.... ');
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
    [h flag] = sample_h(s2,rho,sigh2,h,n);
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

    %% forecast
    %x_t1 = [reshape(shortY(end:-1:end-p+1,:)',1,n*p)];
    parfor iforecast = 1:nforecast
        %[s1, s2] = RandStream.create('mrg32k3a','NumStreams',2); % rng
        x_t1 = [reshape(shortY(end:-1:end-p+1,:)',1,n*p)];
        h_end = h(end);
        for tt = 1:tH
            lam_t1 = 0; mu_t1 = zeros(1,n);
            z_t1 = find(mnrnd(1,[accumarray(z',1);alpha]/(T+alpha))==1);
            if z_t1 > length(lamK)
                lam_t1 = 1./gamrnd(nu0_dp/2,2/S0_dp,1,1);
                mu_t1 = (mu0 + chol(V0/(100^2*lam_t1),'lower')\randn(n,1))';
            else
                lam_t1 = lamK(z_t1);
                mu_t1 = muK(z_t1,:);
            end
            
            h_t1 = rho*h_end + sqrt(sigh2)*randn;
            %mu_t1 + sqrt(lam_t1*exp(h_t1))*(chol(Sig,'lower')*randn(n,1))';
            y_t1 = x_t1*A + mu_t1 + sqrt(lam_t1*exp(h_t1))*(chol(Sig,'lower')*randn(n,1))';
            x_t1 = [y_t1 x_t1(1:end-n)];
            h_end = h_t1;
            store_forecasts(:,tt,iforecast) = y_t1;
        end
    end

    %{
    den1 = zeros(n,1);
    parfor ii=1:n
        pd = fitdist(store_forecasts(:,1,ii),'Kernel','Kernel','normal');
        den1(ii) = pdf(pd,dataT(1,ii));
    end
    
    EY4 = zeros(n,1); den4 = zeros(n,1);
    parfor ii=1:n                
        forecasts4ii = mean(store_forecasts(:,1:4,ii),2);
        Yii4 = mean(dataT(1:4,ii));    
        pd = fitdist(forecasts4ii,'Kernel','Kernel','normal');
        EY4(ii) = mean(forecasts4ii);
        den4(ii) = pdf(pd,Yii4) + 1e-10;
    end 

    EY8 = zeros(n,1); den8 = zeros(n,1);
    parfor ii=1:n                
        forecasts8ii = mean(store_forecasts(:,5:8,ii),2);
        Yii8 = mean(dataT(5:8,ii));    
        pd = fitdist(forecasts8ii,'Kernel','Kernel','normal');
        EY8(ii) = mean(forecasts8ii);
        den8(ii) = pdf(pd,Yii8) + 1e-10;
    end 
    %}
    
    if isim > burnin
        isave = isim - burnin; 
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;
        store_h(isave,:) = h';
        store_lam(isave,:) = lam';
        store_theta(isave,:) = [alpha K rho sigh2];
        store_z(isave,:) = z;
        store_mu(:,:,isave) = mu;

        %tmpyhat1(isave,:) = [squeeze(mean(store_forecasts(:,1,:)))' den1'];
        %tmpyhat4(isave,:) = [EY4' den4'];
        %tmpyhat8(isave,:) = [EY8' den8'];
        %fcstYdraws = zeros(n,tH,nforecast,nsims); 
        fcstYdraws(:,:,:,isave) = store_forecasts;
        
    end
    
    if ( mod(isim, 500) == 0 )
        disp([ num2str(isim) ' loops... ' ])
    end 
    
end

disp( ['MCMC takes '  num2str( etime(clock, start_time) ) 'seconds' ] );
disp(' ');

%{
yhat1 = [dataT(1,:)', mean(tmpyhat1(:,1:n))', log(mean(tmpyhat1(:,n+1:end)))'];
yhat4 = [mean(dataT(1:4,:))', mean(tmpyhat4(:,1:n))', log(mean(tmpyhat4(:,n+1:end)))'];
yhat8 = [mean(dataT(5:8,:))', mean(tmpyhat8(:,1:n))', log(mean(tmpyhat8(:,n+1:end)))'];
store_yhat = [yhat1, yhat4, yhat8];
%}
fcstYdraws = reshape(fcstYdraws, n, tH, nforecast*nsims);

for hh = 1 : tH
    for ii = 1 : n
        fcstCRPS(ii,hh) = crpsDraws(dataT(hh,ii), fcstYdraws(ii,hh,:));
        fcstYhat(ii,hh) = mean(fcstYdraws(ii,hh,:));
    end
end
    
store_yhat = [dataT' fcstYhat fcstCRPS]; % store forecast results.



