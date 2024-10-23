clear; close all; clc;
% 16 vs 99 variable series
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

dpmsv_y1 = NaN(16,24,455);
dpmsv_yhat1 = NaN(16,24,455);
dpmsv_crps1 = NaN(16,24,455);

for tt = 311:765
    str = strcat ('455 dpmsv/',num2str(tt),'_dpmsv99_455.txt');
    dpmsv_data1 = load(str);
    dpmsv_y1(:,:,tt-310) = dpmsv_data1([1 3 5 18 21 29 42 96 76 92 45 59 74 65 66 71],1:24);
    dpmsv_yhat1(:,:,tt-310) = dpmsv_data1([1 3 5 18 21 29 42 96 76 92 45 59 74 65 66 71],25:48);
    dpmsv_crps1(:,:,tt-310) = dpmsv_data1([1 3 5 18 21 29 42 96 76 92 45 59 74 65 66 71],49:72);
end


% check y, dpmsv_y
dy1 = dpmsv_y1 - dpmsv_y;

dy = [sum(sum(sum(dy1,'omitnan'),'omitnan'),'omitnan')];
dy > 1.0e-04

y = dpmsv_y;

dpmsv_mse = NaN(16,24,455);
dpmsv_mse1 = NaN(16,24,455);


for i = 1:455
    dpmsv_mse(:,:,i) = (dpmsv_yhat(:,:,i) - y(:,:,i)).^2;
    dpmsv_mse1(:,:,i) = (dpmsv_yhat1(:,:,i) - y(:,:,i)).^2;
end


%% ======================================================================
% 1985(t=311)-2007(t=586) Great Moderation, 2008(587)-2014(670) GFC,
% 2020:03(733)-2021:02(744) COVID-19
% fcst_T - 310, t311 = 1, t586=276, t587=277, t670=360, t733=423, t744=434

t1 = 277;
t2 = 431; %455-24 = 431 (2020:11)

dpmsv_rmse = sqrt(mean(dpmsv_mse(:,:,t1:t2),3,'omitnan'));
dpmsv_rmse1 = sqrt(mean(dpmsv_mse1(:,:,t1:t2),3,'omitnan'));


dpmsv_acrps = mean(dpmsv_crps(:,:,t1:t2),3,'omitnan');
dpmsv_acrps1 = mean(dpmsv_crps1(:,:,t1:t2),3,'omitnan');



% ratio
ratio_rmse1 = dpmsv_rmse./dpmsv_rmse1;
ratio_rmse1a = mean(ratio_rmse1,2,'omitnan');
ratio_acrps1 = dpmsv_acrps./dpmsv_acrps1;
ratio_acrps1a = mean(ratio_acrps1,2,'omitnan');

results = [ratio_rmse1(:,1) ratio_rmse1(:,3) ratio_rmse1(:,12) ratio_rmse1(:,24)...
    ratio_acrps1(:,1) ratio_acrps1(:,3) ratio_acrps1(:,12) ratio_acrps1(:,24)];

%writematrix(results,'foo/cpp02.csv')












