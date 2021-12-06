clear
clc

cd('F:\My research\case2\JGR\Working\step5_iveg_type')
IDveg = ncread('restart_out_trendy_1_1901_run_C_N_cycle.nc','iveg');
% analysis for vegetation type. Note that only natural vegetation types are included (1-8)
% ID: 1 forest evergreen_needleleaf;
%     2 forest evergreen_broadleaf;
%     3 deciduous_needleleaf;
%     4 deciduous_broadleaf;
%     5 Shrub;
%     6 C3 grassland;
%     7 C4 grassland;
%     8 tundra

% C-only
C_NPP_yrs = [];     % NPP, unit:gC m-2 yr-1, from 1901 to 2013
C_tuaE_yrs = [];    % residence time,      unit: year
C_Xp_yrs = [];      % disequilibrium term, unit: gC m-2
C_X_yrs = [];       % C storage,           unit: gC m-2
C_Xc_yrs = [];      % C storage capacity,  unit: gC m-2

file_NPP = 'F:\My research\case2\working\cable-C\Matrix_out\';
file_tuaE = 'F:\My research\case2\working\cable-C\Transit_Matrix\';   % original unit: day
file_Xp = 'F:\My research\case2\working\cable-C\Transit_Matrix\';
file_X = 'F:\My research\case2\working\cable-C\Transit_Matrix\';
file_Xc = 'F:\My research\case2\working\cable-C\Transit_Matrix\';
ndays = 365;

for i= 1901:2013
    i
    C_NPPinf = xlsread([file_NPP, 'Carbon_1_',num2str(i),'_run_CN_cycle.csv']);   % original unit: gC m-2 yr-1
    C_tuaEinf = xlsread([file_tuaE,'resident_All_',num2str(i),'.csv']);           % original unit: day
    C_Xp9_inf = xlsread([file_Xp, 'C_potencial_',num2str(i),'.csv']);  % the disequilibrium              (unit: gC m-2)
    C_X9_inf = xlsread([file_X, 'X_xx_',num2str(i),'.csv']);           % C storage of 9 C pools          (unit: gC m-2)
    C_Xc9_inf = xlsread([file_Xc,'C_capacity_',num2str(i),'.csv']);    % C storage capacity of 9 C pools (unit: gC m-2)
    
    C_NPP_yrs(:,i-1900) = C_NPPinf(:,1);
    C_tuaE_yrs(:,i-1900) = C_tuaEinf(:,1)./ndays;   % convert unit from day into year
    C_Xp_yrs(:,i-1900) = sum(C_Xp9_inf,2);
    C_X_yrs(:,:,i-1900) = sum(C_X9_inf,2);
    C_Xc_yrs(:,:,i-1900) = sum(C_Xc9_inf,2);
     
end

% calculate the initial values of NPP and residence time
% mean value of the first decade
C_NPP_10bg = nanmean(C_NPP_yrs(:,1:10),2);
C_tuaE_10bg = nanmean(C_tuaE_yrs(:,1:10),2);
C_Xp_10bg = nanmean(C_Xp_yrs(:,1:10),2);
C_X_10bg = nanmean(C_X_yrs(:,1:10),2);
C_Xc_10bg = nanmean(C_Xc_yrs(:,1:10),2);

C_NPP_10ed = nanmean(C_NPP_yrs(:,104:113),2);
C_tuaE_10ed = nanmean(C_tuaE_yrs(:,104:113),2);
C_Xp_10ed = nanmean(C_Xp_yrs(:,104:113),2);
C_X_10ed = nanmean(C_X_yrs(:,104:113),2);
C_Xc_10ed = nanmean(C_Xc_yrs(:,104:113),2);

Net_C_NPP = C_NPP_10ed - C_NPP_10bg;
Net_C_tuaE = C_tuaE_10ed - C_tuaE_10bg;
Net_C_Xp = C_Xp_10ed - C_Xp_10bg;
Net_C_X = C_X_10ed - C_X_10bg;
Net_C_Xc = C_Xc_10ed - C_Xc_10bg;

clearvars -except IDveg C_NPP_yrs C_tuaE_yrs C_Xp_yrs C_X_yrs C_Xc_yrs ...
                  C_NPP_10bg C_tuaE_10bg C_Xp_10bg C_X_10bg C_Xc_10bg ...
                  Net_C_NPP Net_C_tuaE Net_C_Xp Net_C_X Net_C_Xc

% CN
CN_NPP_yrs = [];     % NPP, unit:gC m-2 yr-1, from 1901 to 2013
CN_tuaE_yrs = [];    % residence time,      unit: year
CN_Xp_yrs = [];      % disequilibrium term, unit: gC m-2
CN_X_yrs = [];       % C storage,           unit: gC m-2
CN_Xc_yrs = [];      % C storage capacity,  unit: gC m-2

file_NPP = 'F:\My research\case2\working\cable-CN\Matrix_out\';
file_tuaE = 'F:\My research\case2\working\cable-CN\Transit_Matrix\';   % original unit: day
file_Xp = 'F:\My research\case2\working\cable-CN\Transit_Matrix\';
file_X = 'F:\My research\case2\working\cable-CN\Transit_Matrix\';
file_Xc = 'F:\My research\case2\working\cable-CN\Transit_Matrix\';
ndays = 365;

for i= 1901:2013
    i
    CN_NPPinf = xlsread([file_NPP, 'Carbon_1_',num2str(i),'_run_CN_cycle.csv']);   % original unit: gC m-2 yr-1
    CN_tuaEinf = xlsread([file_tuaE,'resident_All_',num2str(i),'.csv']);           % original unit: day
    CN_Xp9_inf = xlsread([file_Xp, 'C_potencial_',num2str(i),'.csv']);  % the disequilibrium              (unit: gC m-2)
    CN_X9_inf = xlsread([file_X, 'X_xx_',num2str(i),'.csv']);           % C storage of 9 C pools          (unit: gC m-2)
    CN_Xc9_inf = xlsread([file_Xc,'C_capacity_',num2str(i),'.csv']);    % C storage capacity of 9 C pools (unit: gC m-2)
    
    CN_NPP_yrs(:,i-1900) = CN_NPPinf(:,1);
    CN_tuaE_yrs(:,i-1900) = CN_tuaEinf(:,1)./ndays;   % convert unit from day into year
    CN_Xp_yrs(:,i-1900) = sum(CN_Xp9_inf,2);
    CN_X_yrs(:,:,i-1900) = sum(CN_X9_inf,2);
    CN_Xc_yrs(:,:,i-1900) = sum(CN_Xc9_inf,2);
     
end

% calculate the initial values of NPP and residence time
% mean value of the first decade
CN_NPP_10bg = nanmean(CN_NPP_yrs(:,1:10),2);
CN_tuaE_10bg = nanmean(CN_tuaE_yrs(:,1:10),2);
CN_Xp_10bg = nanmean(CN_Xp_yrs(:,1:10),2);
CN_X_10bg = nanmean(CN_X_yrs(:,1:10),2);
CN_Xc_10bg = nanmean(CN_Xc_yrs(:,1:10),2);

CN_NPP_10ed = nanmean(CN_NPP_yrs(:,104:113),2);
CN_tuaE_10ed = nanmean(CN_tuaE_yrs(:,104:113),2);
CN_Xp_10ed = nanmean(CN_Xp_yrs(:,104:113),2);
CN_X_10ed = nanmean(CN_X_yrs(:,104:113),2);
CN_Xc_10ed = nanmean(CN_Xc_yrs(:,104:113),2);

Net_CN_NPP = CN_NPP_10ed - CN_NPP_10bg;
Net_CN_tuaE = CN_tuaE_10ed - CN_tuaE_10bg;
Net_CN_Xp = CN_Xp_10ed - CN_Xp_10bg;
Net_CN_X = CN_X_10ed - CN_X_10bg;
Net_CN_Xc = CN_Xc_10ed - CN_Xc_10bg;

clearvars -except IDveg C_NPP_yrs C_tuaE_yrs C_Xp_yrs C_X_yrs C_Xc_yrs ...
                  C_NPP_10bg C_tuaE_10bg C_Xp_10bg C_X_10bg C_Xc_10bg ...
                  Net_C_NPP Net_C_tuaE Net_C_Xp Net_C_X Net_C_Xc ...
                  CN_NPP_yrs CN_tuaE_yrs CN_Xp_yrs CN_X_yrs CN_Xc_yrs ...
                  CN_NPP_10bg CN_tuaE_10bg CN_Xp_10bg CN_X_10bg CN_Xc_10bg ...
                  Net_CN_NPP Net_CN_tuaE Net_CN_Xp Net_CN_X Net_CN_Xc

% CNP
CNP_NPP_yrs = [];     % NPP, unit:gC m-2 yr-1, from 1901 to 2013
CNP_tuaE_yrs = [];    % residence time,      unit: year
CNP_Xp_yrs = [];      % disequilibrium term, unit: gC m-2
CNP_X_yrs = [];       % C storage,           unit: gC m-2
CNP_Xc_yrs = [];      % C storage capacity,  unit: gC m-2

file_NPP = 'F:\My research\case2\working\cable-CNP\Matrix_out\';
file_tuaE = 'F:\My research\case2\working\cable-CNP\Transit_Matrix\';   % original unit: day
file_Xp = 'F:\My research\case2\working\cable-CNP\Transit_Matrix\';
file_X = 'F:\My research\case2\working\cable-CNP\Transit_Matrix\';
file_Xc = 'F:\My research\case2\working\cable-CNP\Transit_Matrix\';
ndays = 365;

for i= 1901:2013
    i
    CNP_NPPinf = xlsread([file_NPP, 'Carbon_1_',num2str(i),'_run_CN_cycle.csv']);   % original unit: gC m-2 yr-1
    CNP_tuaEinf = xlsread([file_tuaE,'resident_All_',num2str(i),'.csv']);           % original unit: day
    CNP_Xp9_inf = xlsread([file_Xp, 'C_potencial_',num2str(i),'.csv']);  % the disequilibrium              (unit: gC m-2)
    CNP_X9_inf = xlsread([file_X, 'X_xx_',num2str(i),'.csv']);           % C storage of 9 C pools          (unit: gC m-2)
    CNP_Xc9_inf = xlsread([file_Xc,'C_capacity_',num2str(i),'.csv']);    % C storage capacity of 9 C pools (unit: gC m-2)
    
    CNP_NPP_yrs(:,i-1900) = CNP_NPPinf(:,1);
    CNP_tuaE_yrs(:,i-1900) = CNP_tuaEinf(:,1)./ndays;   % convert unit from day into year
    CNP_Xp_yrs(:,i-1900) = sum(CNP_Xp9_inf,2);
    CNP_X_yrs(:,:,i-1900) = sum(CNP_X9_inf,2);
    CNP_Xc_yrs(:,:,i-1900) = sum(CNP_Xc9_inf,2);
     
end

% calculate the initial values of NPP and residence time
% mean value of the first decade
CNP_NPP_10bg = nanmean(CNP_NPP_yrs(:,1:10),2);
CNP_tuaE_10bg = nanmean(CNP_tuaE_yrs(:,1:10),2);
CNP_Xp_10bg = nanmean(CNP_Xp_yrs(:,1:10),2);
CNP_X_10bg = nanmean(CNP_X_yrs(:,1:10),2);
CNP_Xc_10bg = nanmean(CNP_Xc_yrs(:,1:10),2);

CNP_NPP_10ed = nanmean(CNP_NPP_yrs(:,104:113),2);
CNP_tuaE_10ed = nanmean(CNP_tuaE_yrs(:,104:113),2);
CNP_Xp_10ed = nanmean(CNP_Xp_yrs(:,104:113),2);
CNP_X_10ed = nanmean(CNP_X_yrs(:,104:113),2);
CNP_Xc_10ed = nanmean(CNP_Xc_yrs(:,104:113),2);

Net_CNP_NPP = CNP_NPP_10ed - CNP_NPP_10bg;
Net_CNP_tuaE = CNP_tuaE_10ed - CNP_tuaE_10bg;
Net_CNP_Xp = CNP_Xp_10ed - CNP_Xp_10bg;
Net_CNP_X = CNP_X_10ed - CNP_X_10bg;
Net_CNP_Xc = CNP_Xc_10ed - CNP_Xc_10bg;

clearvars -except IDveg C_NPP_yrs C_tuaE_yrs C_Xp_yrs C_X_yrs C_Xc_yrs ...
                  C_NPP_10bg C_tuaE_10bg C_Xp_10bg C_X_10bg C_Xc_10bg ...
                  Net_C_NPP Net_C_tuaE Net_C_Xp Net_C_X Net_C_Xc ...
                  CN_NPP_yrs CN_tuaE_yrs CN_Xp_yrs CN_X_yrs CN_Xc_yrs ...
                  CN_NPP_10bg CN_tuaE_10bg CN_Xp_10bg CN_X_10bg CN_Xc_10bg ...
                  Net_CN_NPP Net_CN_tuaE Net_CN_Xp Net_CN_X Net_CN_Xc ...
                  CNP_NPP_yrs CNP_tuaE_yrs CNP_Xp_yrs CNP_X_yrs CNP_Xc_yrs ...
                  CNP_NPP_10bg CNP_tuaE_10bg CNP_Xp_10bg CNP_X_10bg CNP_Xc_10bg ...
                  Net_CNP_NPP Net_CNP_tuaE Net_CNP_Xp Net_CNP_X Net_CNP_Xc
save ('F:\My research\case2\JGR\1_Figures\Figure4\matData\iveg_reNPPtuaE.mat')
load('F:\My research\case2\JGR\1_Figures\Figure4\matData\iveg_reNPPtuaE.mat')
              
%% classify into differnt veg types
% 1 forest evergreen_needleleaf
index_ENF = find(IDveg==1);

C_ENF_NPP_10bg = C_NPP_10bg(index_ENF);
C_ENF_tuaE_10bg = C_tuaE_10bg(index_ENF);
C_ENF_Xp_10bg = C_Xp_10bg(index_ENF);
C_ENF_X_10bg = C_X_10bg(index_ENF);
C_ENF_Xc_10bg = C_Xc_10bg(index_ENF);

CN_ENF_NPP_10bg = CN_NPP_10bg(index_ENF);
CN_ENF_tuaE_10bg = CN_tuaE_10bg(index_ENF);
CN_ENF_Xp_10bg = CN_Xp_10bg(index_ENF);
CN_ENF_X_10bg = CN_X_10bg(index_ENF);
CN_ENF_Xc_10bg = CN_Xc_10bg(index_ENF);

CNP_ENF_NPP_10bg = CNP_NPP_10bg(index_ENF);
CNP_ENF_tuaE_10bg = CNP_tuaE_10bg(index_ENF);
CNP_ENF_Xp_10bg = CNP_Xp_10bg(index_ENF);
CNP_ENF_X_10bg = CNP_X_10bg(index_ENF);
CNP_ENF_Xc_10bg = CNP_Xc_10bg(index_ENF);

Net_ENF_C_NPP = Net_C_NPP(index_ENF);
Net_ENF_C_tuaE = Net_C_tuaE(index_ENF);
Net_ENF_C_Xp = Net_C_Xp(index_ENF);
Net_ENF_C_X = Net_C_X(index_ENF);
Net_ENF_C_Xc = Net_C_Xc(index_ENF);

Net_ENF_CN_NPP = Net_CN_NPP(index_ENF);
Net_ENF_CN_tuaE = Net_CN_tuaE(index_ENF);
Net_ENF_CN_Xp = Net_CN_Xp(index_ENF);
Net_ENF_CN_X = Net_CN_X(index_ENF);
Net_ENF_CN_Xc = Net_CN_Xc(index_ENF);              
              
Net_ENF_CNP_NPP = Net_CNP_NPP(index_ENF);
Net_ENF_CNP_tuaE = Net_CNP_tuaE(index_ENF);
Net_ENF_CNP_Xp = Net_CNP_Xp(index_ENF);
Net_ENF_CNP_X = Net_CNP_X(index_ENF);
Net_ENF_CNP_Xc = Net_CNP_Xc(index_ENF);               
              
              
% 2 forest evergreen_broadleaf              
index_EBF = find(IDveg==2);  

C_EBF_NPP_10bg = C_NPP_10bg(index_EBF);
C_EBF_tuaE_10bg = C_tuaE_10bg(index_EBF);
C_EBF_Xp_10bg = C_Xp_10bg(index_EBF);
C_EBF_X_10bg = C_X_10bg(index_EBF);
C_EBF_Xc_10bg = C_Xc_10bg(index_EBF);

CN_EBF_NPP_10bg = CN_NPP_10bg(index_EBF);
CN_EBF_tuaE_10bg = CN_tuaE_10bg(index_EBF);
CN_EBF_Xp_10bg = CN_Xp_10bg(index_EBF);
CN_EBF_X_10bg = CN_X_10bg(index_EBF);
CN_EBF_Xc_10bg = CN_Xc_10bg(index_EBF);

CNP_EBF_NPP_10bg = CNP_NPP_10bg(index_EBF);
CNP_EBF_tuaE_10bg = CNP_tuaE_10bg(index_EBF);
CNP_EBF_Xp_10bg = CNP_Xp_10bg(index_EBF);
CNP_EBF_X_10bg = CNP_X_10bg(index_EBF);
CNP_EBF_Xc_10bg = CNP_Xc_10bg(index_EBF);

Net_EBF_C_NPP = Net_C_NPP(index_EBF);
Net_EBF_C_tuaE = Net_C_tuaE(index_EBF);
Net_EBF_C_Xp = Net_C_Xp(index_EBF);
Net_EBF_C_X = Net_C_X(index_EBF);
Net_EBF_C_Xc = Net_C_Xc(index_EBF);

Net_EBF_CN_NPP = Net_CN_NPP(index_EBF);
Net_EBF_CN_tuaE = Net_CN_tuaE(index_EBF);
Net_EBF_CN_Xp = Net_CN_Xp(index_EBF);
Net_EBF_CN_X = Net_CN_X(index_EBF);
Net_EBF_CN_Xc = Net_CN_Xc(index_EBF);              
              
Net_EBF_CNP_NPP = Net_CNP_NPP(index_EBF);
Net_EBF_CNP_tuaE = Net_CNP_tuaE(index_EBF);
Net_EBF_CNP_Xp = Net_CNP_Xp(index_EBF);
Net_EBF_CNP_X = Net_CNP_X(index_EBF);
Net_EBF_CNP_Xc = Net_CNP_Xc(index_EBF);
              
              
% 3 deciduous_needleleaf              
index_DNF = find(IDveg==3);              
              
C_DNF_NPP_10bg = C_NPP_10bg(index_DNF);
C_DNF_tuaE_10bg = C_tuaE_10bg(index_DNF);
C_DNF_Xp_10bg = C_Xp_10bg(index_DNF);
C_DNF_X_10bg = C_X_10bg(index_DNF);
C_DNF_Xc_10bg = C_Xc_10bg(index_DNF);

CN_DNF_NPP_10bg = CN_NPP_10bg(index_DNF);
CN_DNF_tuaE_10bg = CN_tuaE_10bg(index_DNF);
CN_DNF_Xp_10bg = CN_Xp_10bg(index_DNF);
CN_DNF_X_10bg = CN_X_10bg(index_DNF);
CN_DNF_Xc_10bg = CN_Xc_10bg(index_DNF);

CNP_DNF_NPP_10bg = CNP_NPP_10bg(index_DNF);
CNP_DNF_tuaE_10bg = CNP_tuaE_10bg(index_DNF);
CNP_DNF_Xp_10bg = CNP_Xp_10bg(index_DNF);
CNP_DNF_X_10bg = CNP_X_10bg(index_DNF);
CNP_DNF_Xc_10bg = CNP_Xc_10bg(index_DNF);

Net_DNF_C_NPP = Net_C_NPP(index_DNF);
Net_DNF_C_tuaE = Net_C_tuaE(index_DNF);
Net_DNF_C_Xp = Net_C_Xp(index_DNF);
Net_DNF_C_X = Net_C_X(index_DNF);
Net_DNF_C_Xc = Net_C_Xc(index_DNF);

Net_DNF_CN_NPP = Net_CN_NPP(index_DNF);
Net_DNF_CN_tuaE = Net_CN_tuaE(index_DNF);
Net_DNF_CN_Xp = Net_CN_Xp(index_DNF);
Net_DNF_CN_X = Net_CN_X(index_DNF);
Net_DNF_CN_Xc = Net_CN_Xc(index_DNF);              
              
Net_DNF_CNP_NPP = Net_CNP_NPP(index_DNF);
Net_DNF_CNP_tuaE = Net_CNP_tuaE(index_DNF);
Net_DNF_CNP_Xp = Net_CNP_Xp(index_DNF);
Net_DNF_CNP_X = Net_CNP_X(index_DNF);
Net_DNF_CNP_Xc = Net_CNP_Xc(index_DNF);

% 4 deciduous_broadleaf
index_DBF = find(IDveg==4);

C_DBF_NPP_10bg = C_NPP_10bg(index_DBF);
C_DBF_tuaE_10bg = C_tuaE_10bg(index_DBF);
C_DBF_Xp_10bg = C_Xp_10bg(index_DBF);
C_DBF_X_10bg = C_X_10bg(index_DBF);
C_DBF_Xc_10bg = C_Xc_10bg(index_DBF);

CN_DBF_NPP_10bg = CN_NPP_10bg(index_DBF);
CN_DBF_tuaE_10bg = CN_tuaE_10bg(index_DBF);
CN_DBF_Xp_10bg = CN_Xp_10bg(index_DBF);
CN_DBF_X_10bg = CN_X_10bg(index_DBF);
CN_DBF_Xc_10bg = CN_Xc_10bg(index_DBF);

CNP_DBF_NPP_10bg = CNP_NPP_10bg(index_DBF);
CNP_DBF_tuaE_10bg = CNP_tuaE_10bg(index_DBF);
CNP_DBF_Xp_10bg = CNP_Xp_10bg(index_DBF);
CNP_DBF_X_10bg = CNP_X_10bg(index_DBF);
CNP_DBF_Xc_10bg = CNP_Xc_10bg(index_DBF);

Net_DBF_C_NPP = Net_C_NPP(index_DBF);
Net_DBF_C_tuaE = Net_C_tuaE(index_DBF);
Net_DBF_C_Xp = Net_C_Xp(index_DBF);
Net_DBF_C_X = Net_C_X(index_DBF);
Net_DBF_C_Xc = Net_C_Xc(index_DBF);

Net_DBF_CN_NPP = Net_CN_NPP(index_DBF);
Net_DBF_CN_tuaE = Net_CN_tuaE(index_DBF);
Net_DBF_CN_Xp = Net_CN_Xp(index_DBF);
Net_DBF_CN_X = Net_CN_X(index_DBF);
Net_DBF_CN_Xc = Net_CN_Xc(index_DBF);              
              
Net_DBF_CNP_NPP = Net_CNP_NPP(index_DBF);
Net_DBF_CNP_tuaE = Net_CNP_tuaE(index_DBF);
Net_DBF_CNP_Xp = Net_CNP_Xp(index_DBF);
Net_DBF_CNP_X = Net_CNP_X(index_DBF);
Net_DBF_CNP_Xc = Net_CNP_Xc(index_DBF);

% 5 Shrub
index_shrub = find(IDveg==5);

C_shrub_NPP_10bg = C_NPP_10bg(index_shrub);
C_shrub_tuaE_10bg = C_tuaE_10bg(index_shrub);
C_shrub_Xp_10bg = C_Xp_10bg(index_shrub);
C_shrub_X_10bg = C_X_10bg(index_shrub);
C_shrub_Xc_10bg = C_Xc_10bg(index_shrub);

CN_shrub_NPP_10bg = CN_NPP_10bg(index_shrub);
CN_shrub_tuaE_10bg = CN_tuaE_10bg(index_shrub);
CN_shrub_Xp_10bg = CN_Xp_10bg(index_shrub);
CN_shrub_X_10bg = CN_X_10bg(index_shrub);
CN_shrub_Xc_10bg = CN_Xc_10bg(index_shrub);

CNP_shrub_NPP_10bg = CNP_NPP_10bg(index_shrub);
CNP_shrub_tuaE_10bg = CNP_tuaE_10bg(index_shrub);
CNP_shrub_Xp_10bg = CNP_Xp_10bg(index_shrub);
CNP_shrub_X_10bg = CNP_X_10bg(index_shrub);
CNP_shrub_Xc_10bg = CNP_Xc_10bg(index_shrub);

Net_shrub_C_NPP = Net_C_NPP(index_shrub);
Net_shrub_C_tuaE = Net_C_tuaE(index_shrub);
Net_shrub_C_Xp = Net_C_Xp(index_shrub);
Net_shrub_C_X = Net_C_X(index_shrub);
Net_shrub_C_Xc = Net_C_Xc(index_shrub);

Net_shrub_CN_NPP = Net_CN_NPP(index_shrub);
Net_shrub_CN_tuaE = Net_CN_tuaE(index_shrub);
Net_shrub_CN_Xp = Net_CN_Xp(index_shrub);
Net_shrub_CN_X = Net_CN_X(index_shrub);
Net_shrub_CN_Xc = Net_CN_Xc(index_shrub);              
              
Net_shrub_CNP_NPP = Net_CNP_NPP(index_shrub);
Net_shrub_CNP_tuaE = Net_CNP_tuaE(index_shrub);
Net_shrub_CNP_Xp = Net_CNP_Xp(index_shrub);
Net_shrub_CNP_X = Net_CNP_X(index_shrub);
Net_shrub_CNP_Xc = Net_CNP_Xc(index_shrub);

% 6 C3 grassland
index_C3G = find(IDveg==6);

C_C3G_NPP_10bg = C_NPP_10bg(index_C3G);
C_C3G_tuaE_10bg = C_tuaE_10bg(index_C3G);
C_C3G_Xp_10bg = C_Xp_10bg(index_C3G);
C_C3G_X_10bg = C_X_10bg(index_C3G);
C_C3G_Xc_10bg = C_Xc_10bg(index_C3G);

CN_C3G_NPP_10bg = CN_NPP_10bg(index_C3G);
CN_C3G_tuaE_10bg = CN_tuaE_10bg(index_C3G);
CN_C3G_Xp_10bg = CN_Xp_10bg(index_C3G);
CN_C3G_X_10bg = CN_X_10bg(index_C3G);
CN_C3G_Xc_10bg = CN_Xc_10bg(index_C3G);

CNP_C3G_NPP_10bg = CNP_NPP_10bg(index_C3G);
CNP_C3G_tuaE_10bg = CNP_tuaE_10bg(index_C3G);
CNP_C3G_Xp_10bg = CNP_Xp_10bg(index_C3G);
CNP_C3G_X_10bg = CNP_X_10bg(index_C3G);
CNP_C3G_Xc_10bg = CNP_Xc_10bg(index_C3G);

Net_C3G_C_NPP = Net_C_NPP(index_C3G);
Net_C3G_C_tuaE = Net_C_tuaE(index_C3G);
Net_C3G_C_Xp = Net_C_Xp(index_C3G);
Net_C3G_C_X = Net_C_X(index_C3G);
Net_C3G_C_Xc = Net_C_Xc(index_C3G);

Net_C3G_CN_NPP = Net_CN_NPP(index_C3G);
Net_C3G_CN_tuaE = Net_CN_tuaE(index_C3G);
Net_C3G_CN_Xp = Net_CN_Xp(index_C3G);
Net_C3G_CN_X = Net_CN_X(index_C3G);
Net_C3G_CN_Xc = Net_CN_Xc(index_C3G);              
              
Net_C3G_CNP_NPP = Net_CNP_NPP(index_C3G);
Net_C3G_CNP_tuaE = Net_CNP_tuaE(index_C3G);
Net_C3G_CNP_Xp = Net_CNP_Xp(index_C3G);
Net_C3G_CNP_X = Net_CNP_X(index_C3G);
Net_C3G_CNP_Xc = Net_CNP_Xc(index_C3G);

% 7 C4 grassland
index_C4G = find(IDveg==7);

C_C4G_NPP_10bg = C_NPP_10bg(index_C4G);
C_C4G_tuaE_10bg = C_tuaE_10bg(index_C4G);
C_C4G_Xp_10bg = C_Xp_10bg(index_C4G);
C_C4G_X_10bg = C_X_10bg(index_C4G);
C_C4G_Xc_10bg = C_Xc_10bg(index_C4G);

CN_C4G_NPP_10bg = CN_NPP_10bg(index_C4G);
CN_C4G_tuaE_10bg = CN_tuaE_10bg(index_C4G);
CN_C4G_Xp_10bg = CN_Xp_10bg(index_C4G);
CN_C4G_X_10bg = CN_X_10bg(index_C4G);
CN_C4G_Xc_10bg = CN_Xc_10bg(index_C4G);

CNP_C4G_NPP_10bg = CNP_NPP_10bg(index_C4G);
CNP_C4G_tuaE_10bg = CNP_tuaE_10bg(index_C4G);
CNP_C4G_Xp_10bg = CNP_Xp_10bg(index_C4G);
CNP_C4G_X_10bg = CNP_X_10bg(index_C4G);
CNP_C4G_Xc_10bg = CNP_Xc_10bg(index_C4G);

Net_C4G_C_NPP = Net_C_NPP(index_C4G);
Net_C4G_C_tuaE = Net_C_tuaE(index_C4G);
Net_C4G_C_Xp = Net_C_Xp(index_C4G);
Net_C4G_C_X = Net_C_X(index_C4G);
Net_C4G_C_Xc = Net_C_Xc(index_C4G);

Net_C4G_CN_NPP = Net_CN_NPP(index_C4G);
Net_C4G_CN_tuaE = Net_CN_tuaE(index_C4G);
Net_C4G_CN_Xp = Net_CN_Xp(index_C4G);
Net_C4G_CN_X = Net_CN_X(index_C4G);
Net_C4G_CN_Xc = Net_CN_Xc(index_C4G);              
              
Net_C4G_CNP_NPP = Net_CNP_NPP(index_C4G);
Net_C4G_CNP_tuaE = Net_CNP_tuaE(index_C4G);
Net_C4G_CNP_Xp = Net_CNP_Xp(index_C4G);
Net_C4G_CNP_X = Net_CNP_X(index_C4G);
Net_C4G_CNP_Xc = Net_CNP_Xc(index_C4G);

% 8 tundra
index_tdra = find(IDveg==8);

C_tdra_NPP_10bg = C_NPP_10bg(index_tdra);
C_tdra_tuaE_10bg = C_tuaE_10bg(index_tdra);
C_tdra_Xp_10bg = C_Xp_10bg(index_tdra);
C_tdra_X_10bg = C_X_10bg(index_tdra);
C_tdra_Xc_10bg = C_Xc_10bg(index_tdra);

CN_tdra_NPP_10bg = CN_NPP_10bg(index_tdra);
CN_tdra_tuaE_10bg = CN_tuaE_10bg(index_tdra);
CN_tdra_Xp_10bg = CN_Xp_10bg(index_tdra);
CN_tdra_X_10bg = CN_X_10bg(index_tdra);
CN_tdra_Xc_10bg = CN_Xc_10bg(index_tdra);

CNP_tdra_NPP_10bg = CNP_NPP_10bg(index_tdra);
CNP_tdra_tuaE_10bg = CNP_tuaE_10bg(index_tdra);
CNP_tdra_Xp_10bg = CNP_Xp_10bg(index_tdra);
CNP_tdra_X_10bg = CNP_X_10bg(index_tdra);
CNP_tdra_Xc_10bg = CNP_Xc_10bg(index_tdra);

Net_tdra_C_NPP = Net_C_NPP(index_tdra);
Net_tdra_C_tuaE = Net_C_tuaE(index_tdra);
Net_tdra_C_Xp = Net_C_Xp(index_tdra);
Net_tdra_C_X = Net_C_X(index_tdra);
Net_tdra_C_Xc = Net_C_Xc(index_tdra);

Net_tdra_CN_NPP = Net_CN_NPP(index_tdra);
Net_tdra_CN_tuaE = Net_CN_tuaE(index_tdra);
Net_tdra_CN_Xp = Net_CN_Xp(index_tdra);
Net_tdra_CN_X = Net_CN_X(index_tdra);
Net_tdra_CN_Xc = Net_CN_Xc(index_tdra);              
              
Net_tdra_CNP_NPP = Net_CNP_NPP(index_tdra);
Net_tdra_CNP_tuaE = Net_CNP_tuaE(index_tdra);
Net_tdra_CNP_Xp = Net_CNP_Xp(index_tdra);
Net_tdra_CNP_X = Net_CNP_X(index_tdra);
Net_tdra_CNP_Xc = Net_CNP_Xc(index_tdra);


%% Figure4: NPP-driven change ( reNPP*tuaE0) in Xc for different biomes over the period of 1901-2013.
clearvars -except C_ENF_NPP_10bg C_ENF_tuaE_10bg ...
                  Net_ENF_C_NPP Net_ENF_C_tuaE Net_ENF_C_X Net_ENF_C_Xp ...
                  CN_ENF_NPP_10bg CN_ENF_tuaE_10bg ...
                  Net_ENF_CN_NPP Net_ENF_CN_tuaE Net_ENF_CN_X Net_ENF_CN_Xp ...
                  CNP_ENF_NPP_10bg CNP_ENF_tuaE_10bg ...
                  Net_ENF_CNP_NPP Net_ENF_CNP_tuaE Net_ENF_CNP_X Net_ENF_CNP_Xp ...
                  ...
                  C_EBF_NPP_10bg C_EBF_tuaE_10bg ...
                  Net_EBF_C_NPP Net_EBF_C_tuaE Net_EBF_C_X Net_EBF_C_Xp ...
                  CN_EBF_NPP_10bg CN_EBF_tuaE_10bg ...
                  Net_EBF_CN_NPP Net_EBF_CN_tuaE Net_EBF_CN_X Net_EBF_CN_Xp ...
                  CNP_EBF_NPP_10bg CNP_EBF_tuaE_10bg ...
                  Net_EBF_CNP_NPP Net_EBF_CNP_tuaE Net_EBF_CNP_X Net_EBF_CNP_Xp ...
                  ...
                  C_DNF_NPP_10bg C_DNF_tuaE_10bg ...
                  Net_DNF_C_NPP Net_DNF_C_tuaE Net_DNF_C_X Net_DNF_C_Xp ...
                  CN_DNF_NPP_10bg CN_DNF_tuaE_10bg ...
                  Net_DNF_CN_NPP Net_DNF_CN_tuaE Net_DNF_CN_X Net_DNF_CN_Xp ...
                  CNP_DNF_NPP_10bg CNP_DNF_tuaE_10bg ...
                  Net_DNF_CNP_NPP Net_DNF_CNP_tuaE Net_DNF_CNP_X Net_DNF_CNP_Xp ...
                  ...
                  C_DBF_NPP_10bg C_DBF_tuaE_10bg ...
                  Net_DBF_C_NPP Net_DBF_C_tuaE Net_DBF_C_X Net_DBF_C_Xp ...
                  CN_DBF_NPP_10bg CN_DBF_tuaE_10bg ...
                  Net_DBF_CN_NPP Net_DBF_CN_tuaE Net_DBF_CN_X Net_DBF_CN_Xp ...
                  CNP_DBF_NPP_10bg CNP_DBF_tuaE_10bg ...
                  Net_DBF_CNP_NPP Net_DBF_CNP_tuaE Net_DBF_CNP_X Net_DBF_CNP_Xp ...
                  ...
                  C_shrub_NPP_10bg C_shrub_tuaE_10bg ...
                  Net_shrub_C_NPP Net_shrub_C_tuaE Net_shrub_C_X Net_shrub_C_Xp ...
                  CN_shrub_NPP_10bg CN_shrub_tuaE_10bg ...
                  Net_shrub_CN_NPP Net_shrub_CN_tuaE Net_shrub_CN_X Net_shrub_CN_Xp ...
                  CNP_shrub_NPP_10bg CNP_shrub_tuaE_10bg ...
                  Net_shrub_CNP_NPP Net_shrub_CNP_tuaE Net_shrub_CNP_X Net_shrub_CNP_Xp ...
                  ...
                  C_C3G_NPP_10bg C_C3G_tuaE_10bg ...
                  Net_C3G_C_NPP Net_C3G_C_tuaE Net_C3G_C_X Net_C3G_C_Xp ...
                  CN_C3G_NPP_10bg CN_C3G_tuaE_10bg ...
                  Net_C3G_CN_NPP Net_C3G_CN_tuaE Net_C3G_CN_X Net_C3G_CN_Xp ...
                  CNP_C3G_NPP_10bg CNP_C3G_tuaE_10bg ...
                  Net_C3G_CNP_NPP Net_C3G_CNP_tuaE Net_C3G_CNP_X Net_C3G_CNP_Xp ...
                  ...
                  C_C4G_NPP_10bg C_C4G_tuaE_10bg ...
                  Net_C4G_C_NPP Net_C4G_C_tuaE Net_C4G_C_X Net_C4G_C_Xp ...
                  CN_C4G_NPP_10bg CN_C4G_tuaE_10bg ...
                  Net_C4G_CN_NPP Net_C4G_CN_tuaE Net_C4G_CN_X Net_C4G_CN_Xp ...
                  CNP_C4G_NPP_10bg CNP_C4G_tuaE_10bg ...
                  Net_C4G_CNP_NPP Net_C4G_CNP_tuaE Net_C4G_CNP_X Net_C4G_CNP_Xp ...
                  ...
                  C_tdra_NPP_10bg C_tdra_tuaE_10bg ...
                  Net_tdra_C_NPP Net_tdra_C_tuaE Net_tdra_C_X Net_tdra_C_Xp ...
                  CN_tdra_NPP_10bg CN_tdra_tuaE_10bg ...
                  Net_tdra_CN_NPP Net_tdra_CN_tuaE Net_tdra_CN_X Net_tdra_CN_Xp ...
                  CNP_tdra_NPP_10bg CNP_tdra_tuaE_10bg ...
                  Net_tdra_CNP_NPP Net_tdra_CNP_tuaE Net_tdra_CNP_X Net_tdra_CNP_Xp ...
                  Datacell_tuaE0SD_iveg Datacell_NPP0SD_iveg
% NPP0 
Datacell_NPP0M_iveg(1,1) = nanmean(C_ENF_NPP_10bg(:)./1000);   % change unit from gC m-2 yr-1 into KgC m-2 yr-1
Datacell_NPP0M_iveg(1,2) = nanmean(CN_ENF_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(1,3) = nanmean(CNP_ENF_NPP_10bg(:)./1000);

Datacell_NPP0M_iveg(2,1) = nanmean(C_EBF_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(2,2) = nanmean(CN_EBF_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(2,3) = nanmean(CNP_EBF_NPP_10bg(:)./1000);

Datacell_NPP0M_iveg(3,1) = nanmean(C_DNF_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(3,2) = nanmean(CN_DNF_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(3,3) = nanmean(CNP_DNF_NPP_10bg(:)./1000);

Datacell_NPP0M_iveg(4,1) = nanmean(C_DBF_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(4,2) = nanmean(CN_DBF_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(4,3) = nanmean(CNP_DBF_NPP_10bg(:)./1000);

Datacell_NPP0M_iveg(5,1) = nanmean(C_shrub_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(5,2) = nanmean(CN_shrub_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(5,3) = nanmean(CNP_shrub_NPP_10bg(:)./1000);

Datacell_NPP0M_iveg(6,1) = nanmean(C_C3G_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(6,2) = nanmean(CN_C3G_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(6,3) = nanmean(CNP_C3G_NPP_10bg(:)./1000);

Datacell_NPP0M_iveg(7,1) = nanmean(C_C4G_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(7,2) = nanmean(CN_C4G_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(7,3) = nanmean(CNP_C4G_NPP_10bg(:)./1000);

Datacell_NPP0M_iveg(8,1) = nanmean(C_tdra_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(8,2) = nanmean(CN_tdra_NPP_10bg(:)./1000);
Datacell_NPP0M_iveg(8,3) = nanmean(CNP_tdra_NPP_10bg(:)./1000);              
                            
% initial C residence time
Datacell_tuaE0M_iveg(1,1) = nanmean(C_ENF_tuaE_10bg(:));   % unit: year
Datacell_tuaE0M_iveg(1,2) = nanmean(CN_ENF_tuaE_10bg(:));
Datacell_tuaE0M_iveg(1,3) = nanmean(CNP_ENF_tuaE_10bg(:));

Datacell_tuaE0M_iveg(2,1) = nanmean(C_EBF_tuaE_10bg(:));
Datacell_tuaE0M_iveg(2,2) = nanmean(CN_EBF_tuaE_10bg(:));
Datacell_tuaE0M_iveg(2,3) = nanmean(CNP_EBF_tuaE_10bg(:));

Datacell_tuaE0M_iveg(3,1) = nanmean(C_DNF_tuaE_10bg(:));
Datacell_tuaE0M_iveg(3,2) = nanmean(CN_DNF_tuaE_10bg(:));
Datacell_tuaE0M_iveg(3,3) = nanmean(CNP_DNF_tuaE_10bg(:));

Datacell_tuaE0M_iveg(4,1) = nanmean(C_DBF_tuaE_10bg(:));
Datacell_tuaE0M_iveg(4,2) = nanmean(CN_DBF_tuaE_10bg(:));
Datacell_tuaE0M_iveg(4,3) = nanmean(CNP_DBF_tuaE_10bg(:));

Datacell_tuaE0M_iveg(5,1) = nanmean(C_shrub_tuaE_10bg(:));
Datacell_tuaE0M_iveg(5,2) = nanmean(CN_shrub_tuaE_10bg(:));
Datacell_tuaE0M_iveg(5,3) = nanmean(CNP_shrub_tuaE_10bg(:));

Datacell_tuaE0M_iveg(6,1) = nanmean(C_C3G_tuaE_10bg(:));
Datacell_tuaE0M_iveg(6,2) = nanmean(CN_C3G_tuaE_10bg(:));
Datacell_tuaE0M_iveg(6,3) = nanmean(CNP_C3G_tuaE_10bg(:));

Datacell_tuaE0M_iveg(7,1) = nanmean(C_C4G_tuaE_10bg(:));
Datacell_tuaE0M_iveg(7,2) = nanmean(CN_C4G_tuaE_10bg(:));
Datacell_tuaE0M_iveg(7,3) = nanmean(CNP_C4G_tuaE_10bg(:));

Datacell_tuaE0M_iveg(8,1) = nanmean(C_tdra_tuaE_10bg(:));
Datacell_tuaE0M_iveg(8,2) = nanmean(CN_tdra_tuaE_10bg(:));
Datacell_tuaE0M_iveg(8,3) = nanmean(CNP_tdra_tuaE_10bg(:));   

% relative change of Xp (covert unit into KgC m-2)
Datacell_reXp_iveg(1,1) = nanmean(Net_ENF_C_Xp(:)./1000); 
Datacell_reXp_iveg(1,2) = nanmean(Net_ENF_CN_Xp(:)./1000); 
Datacell_reXp_iveg(1,3) = nanmean(Net_ENF_CNP_Xp(:)./1000);

Datacell_reXp_iveg(2,1) = nanmean(Net_EBF_C_Xp(:)./1000);
Datacell_reXp_iveg(2,2) = nanmean(Net_EBF_CN_Xp(:)./1000);
Datacell_reXp_iveg(2,3) = nanmean(Net_EBF_CNP_Xp(:)./1000);

Datacell_reXp_iveg(3,1) = nanmean(Net_DNF_C_Xp(:)./1000);
Datacell_reXp_iveg(3,2) = nanmean(Net_DNF_CN_Xp(:)./1000);
Datacell_reXp_iveg(3,3) = nanmean(Net_DNF_CNP_Xp(:)./1000);

Datacell_reXp_iveg(4,1) = nanmean(Net_DBF_C_Xp(:)./1000);
Datacell_reXp_iveg(4,2) = nanmean(Net_DBF_CN_Xp(:)./1000);
Datacell_reXp_iveg(4,3) = nanmean(Net_DBF_CNP_Xp(:)./1000);

Datacell_reXp_iveg(5,1) = nanmean(Net_shrub_C_Xp(:)./1000);
Datacell_reXp_iveg(5,2) = nanmean(Net_shrub_CN_Xp(:)./1000);
Datacell_reXp_iveg(5,3) = nanmean(Net_shrub_CNP_Xp(:)./1000);

Datacell_reXp_iveg(6,1) = nanmean(Net_C3G_C_Xp(:)./1000);
Datacell_reXp_iveg(6,2) = nanmean(Net_C3G_CN_Xp(:)./1000);
Datacell_reXp_iveg(6,3) = nanmean(Net_C3G_CNP_Xp(:)./1000);

Datacell_reXp_iveg(7,1) = nanmean(Net_C4G_C_Xp(:)./1000);
Datacell_reXp_iveg(7,2) = nanmean(Net_C4G_CN_Xp(:)./1000);
Datacell_reXp_iveg(7,3) = nanmean(Net_C4G_CNP_Xp(:)./1000);

Datacell_reXp_iveg(8,1) = nanmean(Net_tdra_C_Xp(:)./1000);
Datacell_reXp_iveg(8,2) = nanmean(Net_tdra_CN_Xp(:)./1000);
Datacell_reXp_iveg(8,3) = nanmean(Net_tdra_CNP_Xp(:)./1000);
              
% relative change of NPP (unit: gC m-2 yr-1)
Datacell_reNPP_iveg(1,1) = nanmean(Net_ENF_C_NPP(:)./1000); % convert unit from gC m-2 yr-1 into KgC m-2 yr-1
Datacell_reNPP_iveg(1,2) = nanmean(Net_ENF_CN_NPP(:)./1000); 
Datacell_reNPP_iveg(1,3) = nanmean(Net_ENF_CNP_NPP(:)./1000);

Datacell_reNPP_iveg(2,1) = nanmean(Net_EBF_C_NPP(:)./1000);
Datacell_reNPP_iveg(2,2) = nanmean(Net_EBF_CN_NPP(:)./1000);
Datacell_reNPP_iveg(2,3) = nanmean(Net_EBF_CNP_NPP(:)./1000);

Datacell_reNPP_iveg(3,1) = nanmean(Net_DNF_C_NPP(:)./1000);
Datacell_reNPP_iveg(3,2) = nanmean(Net_DNF_CN_NPP(:)./1000);
Datacell_reNPP_iveg(3,3) = nanmean(Net_DNF_CNP_NPP(:)./1000);

Datacell_reNPP_iveg(4,1) = nanmean(Net_DBF_C_NPP(:)./1000);
Datacell_reNPP_iveg(4,2) = nanmean(Net_DBF_CN_NPP(:)./1000);
Datacell_reNPP_iveg(4,3) = nanmean(Net_DBF_CNP_NPP(:)./1000);

Datacell_reNPP_iveg(5,1) = nanmean(Net_shrub_C_NPP(:)./1000);
Datacell_reNPP_iveg(5,2) = nanmean(Net_shrub_CN_NPP(:)./1000);
Datacell_reNPP_iveg(5,3) = nanmean(Net_shrub_CNP_NPP(:)./1000);

Datacell_reNPP_iveg(6,1) = nanmean(Net_C3G_C_NPP(:)./1000);
Datacell_reNPP_iveg(6,2) = nanmean(Net_C3G_CN_NPP(:)./1000);
Datacell_reNPP_iveg(6,3) = nanmean(Net_C3G_CNP_NPP(:)./1000);

Datacell_reNPP_iveg(7,1) = nanmean(Net_C4G_C_NPP(:)./1000);
Datacell_reNPP_iveg(7,2) = nanmean(Net_C4G_CN_NPP(:)./1000);
Datacell_reNPP_iveg(7,3) = nanmean(Net_C4G_CNP_NPP(:)./1000);

Datacell_reNPP_iveg(8,1) = nanmean(Net_tdra_C_NPP(:)./1000);
Datacell_reNPP_iveg(8,2) = nanmean(Net_tdra_CN_NPP(:)./1000);
Datacell_reNPP_iveg(8,3) = nanmean(Net_tdra_CNP_NPP(:)./1000);

%SD
Datacell_reNPP_SD_iveg(1,1) = nanstd(Net_ENF_C_NPP(:)./1000); % convert unit from gC m-2 yr-1 into KgC m-2 yr-1
Datacell_reNPP_SD_iveg(1,2) = nanstd(Net_ENF_CN_NPP(:)./1000); 
Datacell_reNPP_SD_iveg(1,3) = nanstd(Net_ENF_CNP_NPP(:)./1000);

Datacell_reNPP_SD_iveg(2,1) = nanstd(Net_EBF_C_NPP(:)./1000);
Datacell_reNPP_SD_iveg(2,2) = nanstd(Net_EBF_CN_NPP(:)./1000);
Datacell_reNPP_SD_iveg(2,3) = nanstd(Net_EBF_CNP_NPP(:)./1000);

Datacell_reNPP_SD_iveg(3,1) = nanstd(Net_DNF_C_NPP(:)./1000);
Datacell_reNPP_SD_iveg(3,2) = nanstd(Net_DNF_CN_NPP(:)./1000);
Datacell_reNPP_SD_iveg(3,3) = nanstd(Net_DNF_CNP_NPP(:)./1000);

Datacell_reNPP_SD_iveg(4,1) = nanstd(Net_DBF_C_NPP(:)./1000);
Datacell_reNPP_SD_iveg(4,2) = nanstd(Net_DBF_CN_NPP(:)./1000);
Datacell_reNPP_SD_iveg(4,3) = nanstd(Net_DBF_CNP_NPP(:)./1000);

Datacell_reNPP_SD_iveg(5,1) = nanstd(Net_shrub_C_NPP(:)./1000);
Datacell_reNPP_SD_iveg(5,2) = nanstd(Net_shrub_CN_NPP(:)./1000);
Datacell_reNPP_SD_iveg(5,3) = nanstd(Net_shrub_CNP_NPP(:)./1000);

Datacell_reNPP_SD_iveg(6,1) = nanstd(Net_C3G_C_NPP(:)./1000);
Datacell_reNPP_SD_iveg(6,2) = nanstd(Net_C3G_CN_NPP(:)./1000);
Datacell_reNPP_SD_iveg(6,3) = nanstd(Net_C3G_CNP_NPP(:)./1000);

Datacell_reNPP_SD_iveg(7,1) = nanstd(Net_C4G_C_NPP(:)./1000);
Datacell_reNPP_SD_iveg(7,2) = nanstd(Net_C4G_CN_NPP(:)./1000);
Datacell_reNPP_SD_iveg(7,3) = nanstd(Net_C4G_CNP_NPP(:)./1000);

Datacell_reNPP_SD_iveg(8,1) = nanstd(Net_tdra_C_NPP(:)./1000);
Datacell_reNPP_SD_iveg(8,2) = nanstd(Net_tdra_CN_NPP(:)./1000);
Datacell_reNPP_SD_iveg(8,3) = nanstd(Net_tdra_CNP_NPP(:)./1000);

% relative change of tuaE (unit: year)
Datacell_retuaE_iveg(1,1) = nanmean(Net_ENF_C_tuaE(:)); 
Datacell_retuaE_iveg(1,2) = nanmean(Net_ENF_CN_tuaE(:)); 
Datacell_retuaE_iveg(1,3) = nanmean(Net_ENF_CNP_tuaE(:));

Datacell_retuaE_iveg(2,1) = nanmean(Net_EBF_C_tuaE(:));
Datacell_retuaE_iveg(2,2) = nanmean(Net_EBF_CN_tuaE(:));
Datacell_retuaE_iveg(2,3) = nanmean(Net_EBF_CNP_tuaE(:));

Datacell_retuaE_iveg(3,1) = nanmean(Net_DNF_C_tuaE(:));
Datacell_retuaE_iveg(3,2) = nanmean(Net_DNF_CN_tuaE(:));
Datacell_retuaE_iveg(3,3) = nanmean(Net_DNF_CNP_tuaE(:));

Datacell_retuaE_iveg(4,1) = nanmean(Net_DBF_C_tuaE(:));
Datacell_retuaE_iveg(4,2) = nanmean(Net_DBF_CN_tuaE(:));
Datacell_retuaE_iveg(4,3) = nanmean(Net_DBF_CNP_tuaE(:));

Datacell_retuaE_iveg(5,1) = nanmean(Net_shrub_C_tuaE(:));
Datacell_retuaE_iveg(5,2) = nanmean(Net_shrub_CN_tuaE(:));
Datacell_retuaE_iveg(5,3) = nanmean(Net_shrub_CNP_tuaE(:));

Datacell_retuaE_iveg(6,1) = nanmean(Net_C3G_C_tuaE(:));
Datacell_retuaE_iveg(6,2) = nanmean(Net_C3G_CN_tuaE(:));
Datacell_retuaE_iveg(6,3) = nanmean(Net_C3G_CNP_tuaE(:));

Datacell_retuaE_iveg(7,1) = nanmean(Net_C4G_C_tuaE(:));
Datacell_retuaE_iveg(7,2) = nanmean(Net_C4G_CN_tuaE(:));
Datacell_retuaE_iveg(7,3) = nanmean(Net_C4G_CNP_tuaE(:));

Datacell_retuaE_iveg(8,1) = nanmean(Net_tdra_C_tuaE(:));
Datacell_retuaE_iveg(8,2) = nanmean(Net_tdra_CN_tuaE(:));
Datacell_retuaE_iveg(8,3) = nanmean(Net_tdra_CNP_tuaE(:));

% SD
Datacell_retuaE_SD_iveg(1,1) = nanstd(Net_ENF_C_tuaE(:)); 
Datacell_retuaE_SD_iveg(1,2) = nanstd(Net_ENF_CN_tuaE(:)); 
Datacell_retuaE_SD_iveg(1,3) = nanstd(Net_ENF_CNP_tuaE(:));

Datacell_retuaE_SD_iveg(2,1) = nanstd(Net_EBF_C_tuaE(:));
Datacell_retuaE_SD_iveg(2,2) = nanstd(Net_EBF_CN_tuaE(:));
Datacell_retuaE_SD_iveg(2,3) = nanstd(Net_EBF_CNP_tuaE(:));

Datacell_retuaE_SD_iveg(3,1) = nanstd(Net_DNF_C_tuaE(:));
Datacell_retuaE_SD_iveg(3,2) = nanstd(Net_DNF_CN_tuaE(:));
Datacell_retuaE_SD_iveg(3,3) = nanstd(Net_DNF_CNP_tuaE(:));

Datacell_retuaE_SD_iveg(4,1) = nanstd(Net_DBF_C_tuaE(:));
Datacell_retuaE_SD_iveg(4,2) = nanstd(Net_DBF_CN_tuaE(:));
Datacell_retuaE_SD_iveg(4,3) = nanstd(Net_DBF_CNP_tuaE(:));

Datacell_retuaE_SD_iveg(5,1) = nanstd(Net_shrub_C_tuaE(:));
Datacell_retuaE_SD_iveg(5,2) = nanstd(Net_shrub_CN_tuaE(:));
Datacell_retuaE_SD_iveg(5,3) = nanstd(Net_shrub_CNP_tuaE(:));

Datacell_retuaE_SD_iveg(6,1) = nanstd(Net_C3G_C_tuaE(:));
Datacell_retuaE_SD_iveg(6,2) = nanstd(Net_C3G_CN_tuaE(:));
Datacell_retuaE_SD_iveg(6,3) = nanstd(Net_C3G_CNP_tuaE(:));

Datacell_retuaE_SD_iveg(7,1) = nanstd(Net_C4G_C_tuaE(:));
Datacell_retuaE_SD_iveg(7,2) = nanstd(Net_C4G_CN_tuaE(:));
Datacell_retuaE_SD_iveg(7,3) = nanstd(Net_C4G_CNP_tuaE(:));

Datacell_retuaE_SD_iveg(8,1) = nanstd(Net_tdra_C_tuaE(:));
Datacell_retuaE_SD_iveg(8,2) = nanstd(Net_tdra_CN_tuaE(:));
Datacell_retuaE_SD_iveg(8,3) = nanstd(Net_tdra_CNP_tuaE(:));

% relative change of X (covert unit into KgC m-2)
Datacell_reX_iveg(1,1) = nanmean(Net_ENF_C_X(:)./1000); 
Datacell_reX_iveg(1,2) = nanmean(Net_ENF_CN_X(:)./1000); 
Datacell_reX_iveg(1,3) = nanmean(Net_ENF_CNP_X(:)./1000);

Datacell_reX_iveg(2,1) = nanmean(Net_EBF_C_X(:)./1000);
Datacell_reX_iveg(2,2) = nanmean(Net_EBF_CN_X(:)./1000);
Datacell_reX_iveg(2,3) = nanmean(Net_EBF_CNP_X(:)./1000);

Datacell_reX_iveg(3,1) = nanmean(Net_DNF_C_X(:)./1000);
Datacell_reX_iveg(3,2) = nanmean(Net_DNF_CN_X(:)./1000);
Datacell_reX_iveg(3,3) = nanmean(Net_DNF_CNP_X(:)./1000);

Datacell_reX_iveg(4,1) = nanmean(Net_DBF_C_X(:)./1000);
Datacell_reX_iveg(4,2) = nanmean(Net_DBF_CN_X(:)./1000);
Datacell_reX_iveg(4,3) = nanmean(Net_DBF_CNP_X(:)./1000);

Datacell_reX_iveg(5,1) = nanmean(Net_shrub_C_X(:)./1000);
Datacell_reX_iveg(5,2) = nanmean(Net_shrub_CN_X(:)./1000);
Datacell_reX_iveg(5,3) = nanmean(Net_shrub_CNP_X(:)./1000);

Datacell_reX_iveg(6,1) = nanmean(Net_C3G_C_X(:)./1000);
Datacell_reX_iveg(6,2) = nanmean(Net_C3G_CN_X(:)./1000);
Datacell_reX_iveg(6,3) = nanmean(Net_C3G_CNP_X(:)./1000);

Datacell_reX_iveg(7,1) = nanmean(Net_C4G_C_X(:)./1000);
Datacell_reX_iveg(7,2) = nanmean(Net_C4G_CN_X(:)./1000);
Datacell_reX_iveg(7,3) = nanmean(Net_C4G_CNP_X(:)./1000);

Datacell_reX_iveg(8,1) = nanmean(Net_tdra_C_X(:)./1000);
Datacell_reX_iveg(8,2) = nanmean(Net_tdra_CN_X(:)./1000);
Datacell_reX_iveg(8,3) = nanmean(Net_tdra_CNP_X(:)./1000);


save('F:\My research\case2\JGR\1_Figures\Figure4\matData\Figure_iveg_Dataset.mat')
load('F:\My research\case2\JGR\1_Figures\Figure4\matData\Figure_iveg_Dataset.mat')
% Figure4 of reNPP*tuaE0 for different vegetation types
clearvars -except Datacell_NPP0M_iveg Datacell_tuaE0M_iveg ...
                  Datacell_reXp_iveg Datacell_reNPP_iveg Datacell_retuaE_iveg
Data_fct_reNPP = Datacell_reNPP_iveg.*Datacell_tuaE0M_iveg

% panel(a)
figure
set(gcf,'position',[100 100 907 360])
tight_subplot(1,1,[0.03 0.08],[0.13 0.03],[0.075 0.5])
hold on

Bar_fct_reNPP = bar(Data_fct_reNPP)

Bar_fct_reNPP(1).FaceColor = [0.90,0.90,0.90];
Bar_fct_reNPP(2).FaceColor = [0.97,0.85,0.85];
Bar_fct_reNPP(3).FaceColor= [0.85,0.87,0.97];
Bar_fct_reNPP(1).EdgeColor = [0.34 0.33 0.33];
Bar_fct_reNPP(2).EdgeColor = [1 0 0];
Bar_fct_reNPP(3).EdgeColor = [0 0 1];
Bar_fct_reNPP(1).LineWidth = 1.5;
Bar_fct_reNPP(2).LineWidth = 1.5;
Bar_fct_reNPP(3).LineWidth = 1.5;

set(gca,'linewidth',1.2,'box','on')
set(gca,'XLim',[0 9],'YLim',[0 35]);
set(gca,'Fontname','Arial','FontSize',12);
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');

legCNP = legend({'C-only','CN','CNP'});
set(legCNP,'color','none','EdgeColor','none','Fontname','Arial','Fontsize',11)
ylabel('\tau_E_0 X \DeltaNPP (KgC m^-^2) ','Fontname','Arial','FontSize',14)
xticks([1, 2, 3, 4, 5, 6, 7, 8]);
xticklabels({'ENF','EBF','DNF','DBF','shrub','C3G','C4G','tundra'})
text(0.6, 33,'(a)','Fontname','Arial','FontSize',14)

             
% panel(b): reNPP vs tuaE0
sym_iveg = {'*','o','^','s',... % ENF, EBF, DNF, DBF
            'd','v','p','x'}    % Shrub, C3G, C4G, tundra
legIveg_str = {'ENF', 'EBF', 'DNF', 'DBF',...
               'shrub', 'C3G', 'C4G', 'tundra'};
colorCNP = [0.34 0.33 0.33; 1 0 0; 0 0 1];        

tight_subplot(1,1,[0.03 0.1],[0.13 0.03],[0.59 0.08])
hold on

% C-only
for i = 1:8
    plot(Datacell_tuaE0M_iveg(i,1),Datacell_reNPP_iveg(i,1),'Marker',sym_iveg{i},...
     'MarkerEdgeColor', colorCNP(1,:),'LineWidth',1.2,...
     'MarkerFaceColor', colorCNP(1,:),...
     'MarkerSize',11,'LineStyle','none');
         
end

% CN
for i = 1:8
    plot(Datacell_tuaE0M_iveg(i,2),Datacell_reNPP_iveg(i,2),'Marker',sym_iveg{i},...
     'MarkerEdgeColor', colorCNP(2,:),'LineWidth',1.2,...
     'MarkerSize',11,'LineStyle','none');
         
end

% CNP
for i = 1:8
   CNP_iveg(i) = plot(Datacell_tuaE0M_iveg(i,3),Datacell_reNPP_iveg(i,3),'Marker',sym_iveg{i},...
     'MarkerEdgeColor', colorCNP(3,:),'LineWidth',1.2,...
     'MarkerSize',11,'LineStyle','none');
         
end

set(gca,'linewidth',1.2,'box','on')
set(gca,'XLim',[0 600],'YLim',[0 0.3]);
set(gca,'Fontname','Arial','FontSize',12);
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');

X=linspace(0,600,5000)
Y = 30./X;
plot(X,Y,'k:')
text(160,30./160,'30','Fontname','Arial','FontSize',11,...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','center')

Y = 20./X;
plot(X,Y,'k:')
text(400,20./400,'20','Fontname','Arial','FontSize',11,...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','center')

Y = 10./X;
plot(X,Y,'k:')
text(250,10./250,'10','Fontname','Arial','FontSize',11,...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','center')

Y = 5./X;
plot(X,Y,'k:')
text(60,5./60,'5','Fontname','Arial','FontSize',11,...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','center')

Y = 2./X;
plot(X,Y,'k:')
text(80,0.025,'2','Fontname','Arial','FontSize',11,...
    'BackgroundColor',[1 1 1],'HorizontalAlignment','center')

legCNP = legend(CNP_iveg,legIveg_str)
set(legCNP,'color','none','EdgeColor','none','Fontname','Arial','Fontsize',10,'NumColumns',1)
xlabel('\tau_E_0 (year)','Fontname','Arial','FontSize',14)
ylabel('\DeltaNPP (KgC m^-^2 yr^-^1)','Fontname','Arial','FontSize',14)
text(80,0.28,'(b)','Fontname','Arial','FontSize',14,'HorizontalAlignment','center')

% subplot for C-only, CN and CNP
figure
hold on

%err_C_hz = nanstd(Datacell_tuaE0M_iveg(:,1));
%err_C_vt = nanstd(Datacell_reNPP_iveg(:,1));
%errorbar(nanmean(Datacell_tuaE0M_iveg(:,1)),nanmean(Datacell_reNPP_iveg(:,1)),...
%         err_C_hz,'horizontal','.', 'Color','k','LineWidth',1.2)
%errorbar(nanmean(Datacell_tuaE0M_iveg(:,1)),nanmean(Datacell_reNPP_iveg(:,1)),...
%         err_C_vt,'horizontal','.', 'Color','k','LineWidth',1.2)     

plot(nanmean(Datacell_tuaE0M_iveg(:,1)),nanmean(Datacell_reNPP_iveg(:,1)),'Marker','v',...
     'MarkerEdgeColor', colorCNP(1,:),'LineWidth',1.2,...
     'MarkerFaceColor', colorCNP(1,:),...
     'MarkerSize',12,'LineStyle','none');
 
plot(nanmean(Datacell_tuaE0M_iveg(:,2)),nanmean(Datacell_reNPP_iveg(:,2)),'Marker','o',...
     'MarkerEdgeColor', colorCNP(2,:),'LineWidth',1.2,...
     'MarkerFaceColor', colorCNP(2,:),...
     'MarkerSize',12,'LineStyle','none'); 

plot(nanmean(Datacell_tuaE0M_iveg(:,3)),nanmean(Datacell_reNPP_iveg(:,3)),'Marker','o',...
     'MarkerEdgeColor', colorCNP(3,:),'LineWidth',1.2,...
     'MarkerFaceColor', colorCNP(3,:),...
     'MarkerSize',12,'LineStyle','none');  
 
set(gca,'linewidth',1.2,'box','off')
set(gca,'XLim',[145 155],'YLim',[0.04 0.16]);
set(gca,'Fontname','Arial','FontSize',12);
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto'); 

text(153.5,0.155,'C-only','Fontname','Arial','FontSize',12,'HorizontalAlignment','center')
text(147.55,0.09,'CN','Fontname','Arial','FontSize',12,'HorizontalAlignment','center')
text(148.55,0.06,'CNP','Fontname','Arial','FontSize',12,'HorizontalAlignment','center')


%% FigureS: reXp and reNPP*tuaE0
clearvars -except Datacell_NPP0M_iveg Datacell_tuaE0M_iveg ...
                  Datacell_reXp_iveg Datacell_reNPP_iveg Datacell_retuaE_iveg ...
                  Datacell_reXc_iveg Datacell_reX_iveg

% reNPP*tuaE0 + reXp             
Data_fct_reNPP = Datacell_reNPP_iveg.*Datacell_tuaE0M_iveg

figure
hold on

Bar_fct_reNPP = bar(Data_fct_reNPP)

Bar_fct_reNPP(1).FaceColor = [0.90,0.90,0.90];
Bar_fct_reNPP(2).FaceColor = [0.97,0.85,0.85];
Bar_fct_reNPP(3).FaceColor= [0.85,0.87,0.97];

Bar_fct_reNPP(1).EdgeColor = [0.34 0.33 0.33];
Bar_fct_reNPP(2).EdgeColor = [1 0 0];
Bar_fct_reNPP(3).EdgeColor = [0 0 1];

Bar_fct_reNPP(1).LineWidth = 1.5;
Bar_fct_reNPP(2).LineWidth = 1.5;
Bar_fct_reNPP(3).LineWidth = 1.5;


% reXp
Bar_fct_reXp = bar(Datacell_reXp_iveg)

Bar_fct_reXp(1).FaceColor = 'w';
Bar_fct_reXp(2).FaceColor = 'w';
Bar_fct_reXp(3).FaceColor= 'w';

Bar_fct_reXp(1).EdgeColor = [0.34 0.33 0.33];
Bar_fct_reXp(2).EdgeColor = [1 0 0];
Bar_fct_reXp(3).EdgeColor = [0 0 1];

Bar_fct_reXp(1).LineStyle = ':';
Bar_fct_reXp(2).LineStyle = ':';
Bar_fct_reXp(3).LineStyle = ':';

Bar_fct_reXp(1).LineWidth = 1.5;
Bar_fct_reXp(2).LineWidth = 1.5;
Bar_fct_reXp(3).LineWidth = 1.5;

set(gca,'linewidth',1.2,'box','on')
set(gca,'XLim',[0 9],'YLim',[0 35]);
set(gca,'Fontname','Arial','FontSize',12);
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');

legCNP = legend({'C-only','CN','CNP','Xp'});
set(legCNP,'color','none','EdgeColor','none','Fontname','Arial','Fontsize',11)
ylabel('Ecosystem C change (KgC m^-^2) ','Fontname','Arial','FontSize',14)
xticks([1, 2, 3, 4, 5, 6, 7, 8]);
xticklabels({'ENF','EBF','DNF','DBF','shrub','C3G','C4G','tundra'})
%text(0.6, 33,'(c)','Fontname','Arial','FontSize',14)



