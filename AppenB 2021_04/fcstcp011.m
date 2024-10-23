clear; close all; clc;
% ccmm16 2021-04 T=745 ccmm series
% 1975(t=191)-1984(t=310) Great Inflation,
% 1985(t=311)-2007(t=586) Great Moderation,
% 2008(587)-2014(670) GFC,
% 2020:03(733)-2021:02(744) COVID-19,
% fcst_T - 190, t191 = 1, t586=396, t587=397, t670=480, t733=543, t744=554

dpmsv_y = NaN(16,24,555);
dpmsv_yhat = NaN(16,24,555);
dpmsv_crps = NaN(16,24,555);

for tt = 191:745
    str = strcat ('dpmsv 4p/',num2str(tt),'_ccmm16_A4.txt');
    dpmsv_data = load(str);
    dpmsv_y(:,:,tt-190) = dpmsv_data(:,1:24);
    dpmsv_yhat(:,:,tt-190) = dpmsv_data(:,25:48);
    dpmsv_crps(:,:,tt-190) = dpmsv_data(:,49:72);
end

SVo_y1 = NaN(16,24,555);
SVo_yhat1 = NaN(16,24,555);
SVo_crps1 = NaN(16,24,555);

SVo_data11 = load("fredMD16-2021-04-censoredYields-SVobarmax20-1975-p12.mat");
SVo_data21 = load("fredMD16-2021-04-censoredYields-SVobarmax20-p12.mat");

for i = 1:120
    SVo_y1(:,:,i) = SVo_data11.fcstYrealized(:,:,i);
    SVo_yhat1(:,:,i) = SVo_data11.fcstYhat(:,:,i);
    SVo_crps1(:,:,i) = SVo_data11.fcstCRPS(:,:,i);
end

for i = 121:555
    SVo_y1(:,:,i) = SVo_data21.fcstYrealized(:,:,i-120);
    SVo_yhat1(:,:,i) = SVo_data21.fcstYhat(:,:,i-120);
    SVo_crps1(:,:,i) = SVo_data21.fcstCRPS(:,:,i-120);
end

SVO_y1 = NaN(16,24,555);
SVO_yhat1 = NaN(16,24,555);
SVO_crps1 = NaN(16,24,555);

SVO_data1 = load("fredMD16-2021-04-censoredYields-SVOmax20-p12.mat");

for i = 1:555
    SVO_y1(:,:,i) = SVO_data1.fcstYrealized(:,:,i);
    SVO_yhat1(:,:,i) = SVO_data1.fcstYhat(:,:,i);
    SVO_crps1(:,:,i) = SVO_data1.fcstCRPS(:,:,i);
end

SVo_y = NaN(16,24,555);
SVo_yhat = NaN(16,24,555);
SVo_crps = NaN(16,24,555);

SVo_data1 = load("fredMD16-2021-04-NOshadowrate-SVobarmax20-1975-p12.mat");
SVo_data2 = load("fredMD16-2021-04-NOshadowrate-SVobarmax20-p12.mat");

for i = 1:120
    SVo_y(:,:,i) = SVo_data1.fcstYrealized(:,:,i);
    SVo_yhat(:,:,i) = SVo_data1.fcstYhat(:,:,i);
    SVo_crps(:,:,i) = SVo_data1.fcstCRPS(:,:,i);
end

for i = 121:555
    SVo_y(:,:,i) = SVo_data2.fcstYrealized(:,:,i-120);
    SVo_yhat(:,:,i) = SVo_data2.fcstYhat(:,:,i-120);
    SVo_crps(:,:,i) = SVo_data2.fcstCRPS(:,:,i-120);
end

SVO_y = NaN(16,24,555);
SVO_yhat = NaN(16,24,555);
SVO_crps = NaN(16,24,555);

SVO_data = load("fredMD16-2021-04-NOshadowrate-SVOmax20-p12.mat");

for i = 1:555
    SVO_y(:,:,i) = SVO_data.fcstYrealized(:,:,i);
    SVO_yhat(:,:,i) = SVO_data.fcstYhat(:,:,i);
    SVO_crps(:,:,i) = SVO_data.fcstCRPS(:,:,i);
end

C_y = NaN(16,24,555);
C_yhat = NaN(16,24,555);
C_crps = NaN(16,24,555);

C_data = load("fredMD16-2021-04-NOshadowrate-CONST-p12.mat");

for i = 1:555
    C_y(:,:,i) = C_data.fcstYrealized(:,:,i);
    C_yhat(:,:,i) = C_data.fcstYhat(:,:,i);
    C_crps(:,:,i) = C_data.fcstCRPS(:,:,i);
end

SV_y = NaN(16,24,555);
SV_yhat = NaN(16,24,555);
SV_crps = NaN(16,24,555);

SV_data = load("fredMD16-2021-04-NOshadowrate-SV-p12.mat");

for i = 1:555
    SV_y(:,:,i) = SV_data.fcstYrealized(:,:,i);
    SV_yhat(:,:,i) = SV_data.fcstYhat(:,:,i);
    SV_crps(:,:,i) = SV_data.fcstCRPS(:,:,i);
end

% check y, dpmsv_y
dy1 = SVo_y - dpmsv_y;
dy2 = SVO_y - dpmsv_y;
dy3 = C_y - dpmsv_y;
dy4 = SV_y - dpmsv_y;
dy5 = SVo_y1 - dpmsv_y;
dy6 = SVO_y1 - dpmsv_y;
dy = [sum(sum(sum(dy1,'omitnan'),'omitnan'),'omitnan');...
        sum(sum(sum(dy2,'omitnan'),'omitnan'),'omitnan');...
        sum(sum(sum(dy3,'omitnan'),'omitnan'),'omitnan');...
        sum(sum(sum(dy4,'omitnan'),'omitnan'),'omitnan');...
        sum(sum(sum(dy5,'omitnan'),'omitnan'),'omitnan');...
        sum(sum(sum(dy6,'omitnan'),'omitnan'),'omitnan')];
dy > 1.0e-04 % check the stored true y in fcst result files are the same.

y = dpmsv_y;

dpmsv_mse = NaN(16,24,555);
SVo_mse = NaN(16,24,555);
SVO_mse = NaN(16,24,555);
C_mse = NaN(16,24,555);
SV_mse = NaN(16,24,555);
SVo_mse1 = NaN(16,24,555);
SVO_mse1 = NaN(16,24,555);

for i = 1:555
    dpmsv_mse(:,:,i) = (dpmsv_yhat(:,:,i) - y(:,:,i)).^2;
    SVo_mse(:,:,i) = (SVo_yhat(:,:,i) - y(:,:,i)).^2;
    SVO_mse(:,:,i) = (SVO_yhat(:,:,i) - y(:,:,i)).^2;
    C_mse(:,:,i) = (C_yhat(:,:,i) - y(:,:,i)).^2;
    SV_mse(:,:,i) = (SV_yhat(:,:,i) - y(:,:,i)).^2;
    SVo_mse1(:,:,i) = (SVo_yhat1(:,:,i) - y(:,:,i)).^2;
    SVO_mse1(:,:,i) = (SVO_yhat1(:,:,i) - y(:,:,i)).^2;
end

%{
dpmsv_rmse = NaN(16,24);
SVo_rmse = NaN(16,24);
SVO_rmse = NaN(16,24);
C_rmse = NaN(16,24);
SV_rmse = NaN(16,24);
SVo_rmse1 = NaN(16,24);
SVO_rmse1 = NaN(16,24);

dpmsv_acrps = NaN(16,24);
SVo_acrps = NaN(16,24);
SVO_acrps = NaN(16,24);
C_acrps = NaN(16,24);
SV_acrps = NaN(16,24);
SVo_acrps1 = NaN(16,24);
SVO_acrps1 = NaN(16,24);
%}

%% =================================================================
% 1975(t=191)-1984(t=310) Great Inflation
% 1985(t=311)-2007(t=586) Great Moderation, 2008(587)-2014(670) GFC,
% 2020:03(733)-2021:02(744) COVID-19
% fcst_T - 190, t191=1, t310=120, t311=121, t586=396, 
% t587=397, t670=480, t733=543, t744=554

t1 = 1;
t2 = 120; % 555-24 = 531

dpmsv_rmse = sqrt(mean(dpmsv_mse(:,:,t1:t2),3,'omitnan'));
SVo_rmse = sqrt(mean(SVo_mse(:,:,t1:t2),3,'omitnan'));
SVO_rmse = sqrt(mean(SVO_mse(:,:,t1:t2),3,'omitnan'));
C_rmse = sqrt(mean(C_mse(:,:,t1:t2),3,'omitnan'));
SV_rmse = sqrt(mean(SV_mse(:,:,t1:t2),3,'omitnan'));
SVo_rmse1 = sqrt(mean(SVo_mse1(:,:,t1:t2),3,'omitnan'));
SVO_rmse1 = sqrt(mean(SVO_mse1(:,:,t1:t2),3,'omitnan'));

dpmsv_acrps = mean(dpmsv_crps(:,:,t1:t2),3,'omitnan');
SVo_acrps = mean(SVo_crps(:,:,t1:t2),3,'omitnan');
SVO_acrps = mean(SVO_crps(:,:,t1:t2),3,'omitnan');
C_acrps = mean(C_crps(:,:,t1:t2),3,'omitnan');
SV_acrps = mean(SV_crps(:,:,t1:t2),3,'omitnan');
SVo_acrps1 = mean(SVo_crps1(:,:,t1:t2),3,'omitnan');
SVO_acrps1 = mean(SVO_crps1(:,:,t1:t2),3,'omitnan');

% ratio
ratio_rmse1 = dpmsv_rmse./SVo_rmse;
%ratio_rmse1a = mean(ratio_rmse1,2,'omitnan');
ratio_acrps1 = dpmsv_acrps./SVo_acrps;
%ratio_acrps1a = mean(ratio_acrps1,2,'omitnan');

ratio_rmse2 = dpmsv_rmse./SVO_rmse;
%ratio_rmse2a = mean(ratio_rmse2,2,'omitnan');
ratio_acrps2 = dpmsv_acrps./SVO_acrps;
%ratio_acrps2a = mean(ratio_acrps2,2,'omitnan');

ratio_rmse3 = dpmsv_rmse./C_rmse;
%ratio_rmse3a = mean(ratio_rmse3,2,'omitnan');
ratio_acrps3 = dpmsv_acrps./C_acrps;
%ratio_acrps3a = mean(ratio_acrps3,2,'omitnan');

ratio_rmse4 = dpmsv_rmse./SV_rmse;
%ratio_rmse4a = mean(ratio_rmse4,2,'omitnan');
ratio_acrps4 = dpmsv_acrps./SV_acrps;
%ratio_acrps4a = mean(ratio_acrps4,2,'omitnan');

ratio_rmse5 = dpmsv_rmse./SVo_rmse1;
%ratio_rmse5a = mean(ratio_rmse5,2,'omitnan');
ratio_acrps5 = dpmsv_acrps./SVo_acrps1;
%ratio_acrps5a = mean(ratio_acrps5,2,'omitnan');

ratio_rmse6 = dpmsv_rmse./SVO_rmse1;
%ratio_rmse6a = mean(ratio_rmse6,2,'omitnan');
ratio_acrps6 = dpmsv_acrps./SVO_acrps1;
%ratio_acrps6a = mean(ratio_acrps6,2,'omitnan');

results = [ratio_rmse1(:,1) ratio_rmse1(:,3) ratio_rmse1(:,12) ratio_rmse1(:,24)...
    ratio_acrps1(:,1) ratio_acrps1(:,3) ratio_acrps1(:,12) ratio_acrps1(:,24)...
    ratio_rmse2(:,1) ratio_rmse2(:,3) ratio_rmse2(:,12) ratio_rmse2(:,24)...
    ratio_acrps2(:,1) ratio_acrps2(:,3) ratio_acrps2(:,12) ratio_acrps2(:,24)...
    ratio_rmse3(:,1) ratio_rmse3(:,3) ratio_rmse3(:,12) ratio_rmse3(:,24)...
    ratio_acrps3(:,1) ratio_acrps3(:,3) ratio_acrps3(:,12) ratio_acrps3(:,24)...
    ratio_rmse4(:,1) ratio_rmse4(:,3) ratio_rmse4(:,12) ratio_rmse4(:,24)...
    ratio_acrps4(:,1) ratio_acrps4(:,3) ratio_acrps4(:,12) ratio_acrps4(:,24)...
    ratio_rmse5(:,1) ratio_rmse5(:,3) ratio_rmse5(:,12) ratio_rmse5(:,24)...
    ratio_acrps5(:,1) ratio_acrps5(:,3) ratio_acrps5(:,12) ratio_acrps5(:,24)...
    ratio_rmse6(:,1) ratio_rmse6(:,3) ratio_rmse6(:,12) ratio_rmse6(:,24)...
    ratio_acrps6(:,1) ratio_acrps6(:,3) ratio_acrps6(:,12) ratio_acrps6(:,24)];

writematrix(results,'foo/ccmm_cp01.csv')












