clear; close all; clc;
% ccmm16 2022-12 T=765 return series
% 1985(t=311)-2007(t=586) Great Moderation, 2008(587)-2014(670) GFC,
% 2020:03(733)-2021:02(744) COVID-19
% fcst_T - 310, t311 = 1, t586=276, t587=277, t670=360, t733=423, t744=434

dpmsv_y = NaN(16,24,455);
dpmsv_yhat = NaN(16,24,455);
dpmsv_crps = NaN(16,24,455);

for tt = 311:765
    str = strcat ('455 dpmsv 16/',num2str(tt),'_dpmsv16_455.txt');
    dpmsv_data = load(str);
    dpmsv_y(:,:,tt-310) = dpmsv_data(:,1:24);
    dpmsv_yhat(:,:,tt-310) = dpmsv_data(:,25:48);
    dpmsv_crps(:,:,tt-310) = dpmsv_data(:,49:72);
end

SVo_y = NaN(16,24,455);
SVo_yhat = NaN(16,24,455);
SVo_crps = NaN(16,24,455);

% 2012:01 t635 fcst_T=325
SVo_data1 = load("fredMD16-2022-12-r-NOshadowrate-SVobarmax20-p12.mat");
SVo_data2 = load("fredMD16-2022-12-r-NOshadowrate-SVobarmax20-2008-p12.mat");
SVo_data3 = load("fredMD16-2022-12-r-NOshadowrate-SVobarmax20-1985-p12.mat");

for i = 1:131
    SVo_y(:,:,i+324) = SVo_data1.fcstYrealized(:,:,i);
    SVo_yhat(:,:,i+324) = SVo_data1.fcstYhat(:,:,i);
    SVo_crps(:,:,i+324) = SVo_data1.fcstCRPS(:,:,i);
end

for i = 1:48
    SVo_y(:,:,i+276) = SVo_data2.fcstYrealized(:,:,i);
    SVo_yhat(:,:,i+276) = SVo_data2.fcstYhat(:,:,i);
    SVo_crps(:,:,i+276) = SVo_data2.fcstCRPS(:,:,i);
end

for i = 1:276
    SVo_y(:,:,i) = SVo_data3.fcstYrealized(:,:,i);
    SVo_yhat(:,:,i) = SVo_data3.fcstYhat(:,:,i);
    SVo_crps(:,:,i) = SVo_data3.fcstCRPS(:,:,i);
end

SVO_y = NaN(16,24,455);
SVO_yhat = NaN(16,24,455);
SVO_crps = NaN(16,24,455);

SVO_data = load("fredMD16-2022-12-r-NOshadowrate-SVOmax20-p12.mat");

for i = 1:455
    SVO_y(:,:,i) = SVO_data.fcstYrealized(:,:,i);
    SVO_yhat(:,:,i) = SVO_data.fcstYhat(:,:,i);
    SVO_crps(:,:,i) = SVO_data.fcstCRPS(:,:,i);
end

C_y = NaN(16,24,455);
C_yhat = NaN(16,24,455);
C_crps = NaN(16,24,455);

C_data = load("fredMD16-2022-12-r-NOshadowrate-CONST-p12.mat");

for i = 1:455
    C_y(:,:,i) = C_data.fcstYrealized(:,:,i);
    C_yhat(:,:,i) = C_data.fcstYhat(:,:,i);
    C_crps(:,:,i) = C_data.fcstCRPS(:,:,i);
end

SV_y = NaN(16,24,455);
SV_yhat = NaN(16,24,455);
SV_crps = NaN(16,24,455);

% 2008:01 t635 fcst_T=325
SV_data1 = load("fredMD16-2022-12-r-NOshadowrate-SV-1985-p12.mat");
SV_data2 = load("fredMD16-2022-12-r-NOshadowrate-SV-p12.mat");


for i = 1:276
    SV_y(:,:,i) = SV_data1.fcstYrealized(:,:,i);
    SV_yhat(:,:,i) = SV_data1.fcstYhat(:,:,i);
    SVo_crps(:,:,i) = SV_data1.fcstCRPS(:,:,i);
end

for i = 277:455
    SV_y(:,:,i) = SV_data2.fcstYrealized(:,:,i-276);
    SV_yhat(:,:,i) = SV_data2.fcstYhat(:,:,i-276);
    SV_crps(:,:,i) = SV_data2.fcstCRPS(:,:,i-276);
end

% check y, dpmsv_y
dy1 = SVo_y - dpmsv_y;
dy2 = SVO_y - dpmsv_y;
dy3 = C_y - dpmsv_y;
dy4 = SV_y - dpmsv_y;
dy = [sum(sum(sum(dy1,'omitnan'),'omitnan'),'omitnan');...
        sum(sum(sum(dy2,'omitnan'),'omitnan'),'omitnan');...
        sum(sum(sum(dy3,'omitnan'),'omitnan'),'omitnan');...
        sum(sum(sum(dy4,'omitnan'),'omitnan'),'omitnan')];
dy > 1.0e-04

y = dpmsv_y;

dpmsv_mse = NaN(16,24,455);
SVo_mse = NaN(16,24,455);
SVO_mse = NaN(16,24,455);
C_mse = NaN(16,24,455);
SV_mse = NaN(16,24,455);

for i = 1:455
    dpmsv_mse(:,:,i) = (dpmsv_yhat(:,:,i) - y(:,:,i)).^2;
    SVo_mse(:,:,i) = (SVo_yhat(:,:,i) - y(:,:,i)).^2;
    SVO_mse(:,:,i) = (SVO_yhat(:,:,i) - y(:,:,i)).^2;
    C_mse(:,:,i) = (C_yhat(:,:,i) - y(:,:,i)).^2;
    SV_mse(:,:,i) = (SV_yhat(:,:,i) - y(:,:,i)).^2;
end

%{
dpmsv_rmse = NaN(16,24);
SVo_rmse = NaN(16,24);
SVO_rmse = NaN(16,24);
C_rmse = NaN(16,24);

dpmsv_acrps = NaN(16,24);
SVo_acrps = NaN(16,24);
SVO_acrps = NaN(16,24);
C_acrps = NaN(16,24);
%}

%% ======================================================================
% 1985(t=311)-2007(t=586) Great Moderation, 2008(587)-2014(670) GFC,
% 2020:03(733)-2021:02(744) COVID-19
% fcst_T - 310, t311 = 1, t586=276, t587=277, t670=360, t733=423, t744=434

t1 = 277;
t2 = 434; %455-24 = 431 (2020:11)

dpmsv_rmse = sqrt(mean(dpmsv_mse(:,:,t1:t2),3,'omitnan'));
SVo_rmse = sqrt(mean(SVo_mse(:,:,t1:t2),3,'omitnan'));
SVO_rmse = sqrt(mean(SVO_mse(:,:,t1:t2),3,'omitnan'));
C_rmse = sqrt(mean(C_mse(:,:,t1:t2),3,'omitnan'));
SV_rmse = sqrt(mean(SV_mse(:,:,t1:t2),3,'omitnan'));

dpmsv_acrps = mean(dpmsv_crps(:,:,t1:t2),3,'omitnan');
SVo_acrps = mean(SVo_crps(:,:,t1:t2),3,'omitnan');
SVO_acrps = mean(SVO_crps(:,:,t1:t2),3,'omitnan');
C_acrps = mean(C_crps(:,:,t1:t2),3,'omitnan');
SV_acrps = mean(SV_crps(:,:,t1:t2),3,'omitnan');


% ratio
ratio_rmse1 = dpmsv_rmse./SVo_rmse;
ratio_rmse1a = mean(ratio_rmse1,2,'omitnan');
ratio_acrps1 = dpmsv_acrps./SVo_acrps;
ratio_acrps1a = mean(ratio_acrps1,2,'omitnan');

ratio_rmse2 = dpmsv_rmse./SVO_rmse;
ratio_rmse2a = mean(ratio_rmse2,2,'omitnan');
ratio_acrps2 = dpmsv_acrps./SVO_acrps;
ratio_acrps2a = mean(ratio_acrps2,2,'omitnan');

ratio_rmse3 = dpmsv_rmse./C_rmse;
ratio_rmse3a = mean(ratio_rmse3,2,'omitnan');
ratio_acrps3 = dpmsv_acrps./C_acrps;
ratio_acrps3a = mean(ratio_acrps3,2,'omitnan');

ratio_rmse4 = dpmsv_rmse./SV_rmse;
ratio_rmse4a = mean(ratio_rmse4,2,'omitnan');
ratio_acrps4 = dpmsv_acrps./SV_acrps;
ratio_acrps4a = mean(ratio_acrps4,2,'omitnan');

results = [ratio_rmse1(:,1) ratio_rmse1(:,3) ratio_rmse1(:,12) ratio_rmse1(:,24)...
    ratio_acrps1(:,1) ratio_acrps1(:,3) ratio_acrps1(:,12) ratio_acrps1(:,24)...
    ratio_rmse2(:,1) ratio_rmse2(:,3) ratio_rmse2(:,12) ratio_rmse2(:,24)...
    ratio_acrps2(:,1) ratio_acrps2(:,3) ratio_acrps2(:,12) ratio_acrps2(:,24)...
    ratio_rmse3(:,1) ratio_rmse3(:,3) ratio_rmse3(:,12) ratio_rmse3(:,24)...
    ratio_acrps3(:,1) ratio_acrps3(:,3) ratio_acrps3(:,12) ratio_acrps3(:,24)...
    ratio_rmse4(:,1) ratio_rmse4(:,3) ratio_rmse4(:,12) ratio_rmse4(:,24)...
    ratio_acrps4(:,1) ratio_acrps4(:,3) ratio_acrps4(:,12) ratio_acrps4(:,24)];

%writematrix(results,'foo/cp03.csv')










