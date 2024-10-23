% This file is used to automatically generate individual files running in the cluster
% '2022-12-ccmm16-r.txt', 3/1959 to 11/2022. T = 765
% T1 = 311, 1/1985
% '2022-12-r.txt' full 99 variables
clc;
clear all; 
T1 = 311;
for tt = T1:765
    str = ['forecast_file',num2str(tt),'.m'];
    fd = fopen(str,'a+');
    fprintf(fd,'clc;\r\n');
    fprintf(fd,'clear all;\r\n');
    fprintf(fd,'disp(''file%d is currently running'');\r\n',tt);
    fprintf(fd,'totaltime = tic;\r\n');
    fprintf(fd,'data = readmatrix(''2022-12-r.txt'');\r\n');
    fprintf(fd,'data = [data;NaN(24,size(data,2))];\r\n');
    fprintf(fd,'dataT = data(%d+1:%d+24,:);\r\n',tt,tt); 
    fprintf(fd,'data0 = data(1:%d,:);\r\n',tt);  

    fprintf(fd,'BVAR_DP_CSV_fred_008;\r\n');           
    
    fprintf(fd,'dlmwrite([''%d_dpmsv99_455.txt''],[store_yhat],''precision'',50);\r\n',tt);
    fprintf(fd,'totaltimespend = toc(totaltime)');
    fclose(fd);
end

