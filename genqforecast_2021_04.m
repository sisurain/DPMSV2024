% This file is used to automatically generate individual files running in the cluster
% 'fredMD16-2021-04-nvn.csv', 3/1959 to 3/2021. T = 745
% T1 = 191, 1/1975
clc;
clear all; 
T1 = 191;
for tt = T1:745
    str = ['forecast_file',num2str(tt),'.m'];
    fd = fopen(str,'a+');
    fprintf(fd,'clc;\r\n');
    fprintf(fd,'clear all;\r\n');
    fprintf(fd,'disp(''file%d is currently running'');\r\n',tt);
    fprintf(fd,'totaltime = tic;\r\n');
    fprintf(fd,'data = readmatrix(''fredMD16-2021-04-nvn.csv'');\r\n');
    fprintf(fd,'data = [data;NaN(24,size(data,2))];\r\n');
    fprintf(fd,'dataT = data(%d+1:%d+24,:);\r\n',tt,tt); 
    fprintf(fd,'data0 = data(1:%d,:);\r\n',tt);  
    %fprintf(fd,'data0 = data0(:,[2 16 19 31 32 56 77 90 117 140 144 167 178 180 192 219]);\r\n');
    %fprintf(fd,'dataT = dataT(:,[2 16 19 31 32 56 77 90 117 140 144 167 178 180 192 219]);\r\n');
    fprintf(fd,'BVAR_DP_CSV_fred_008;\r\n');           
    
    fprintf(fd,'dlmwrite([''%d_ccmm16_A4.txt''],[store_yhat],''precision'',50);\r\n',tt);
    fprintf(fd,'totaltimespend = toc(totaltime)');
    fclose(fd);
end

