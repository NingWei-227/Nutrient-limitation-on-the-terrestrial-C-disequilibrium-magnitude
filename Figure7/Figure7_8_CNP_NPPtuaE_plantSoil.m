clear
clc

% 1. calculate C input and residence time for both plant and soil

% read mask and cellarea data used to convert model outputs into global map
[num,text,raw]=xlsread('G:\My research\case2\working\ALL\codes\yanmo0.5.csv');
for rowID=1:360
    for colID =1:720
        if raw{rowID,colID}=='NA'
            raw{rowID,colID}=-9999;
        end  
    end
end
mask = cell2mat(raw);
mask(mask==-9999) =  nan;
cellarea = xlsread('G:\My research\case2\working\ALL\codes\area0.5.csv');   % unit: km2 = 10^6 m2
clearvars -except mask cellarea

%% C-only
% plant and soil
Cp_NPP_maps = [];
Cp_tuaE_maps = [];
Cp_Xc_maps = [];
Cp_Xp_maps = [];

Cs_inC_maps = [];
Cs_tuaE_maps = [];
Cs_Xc_maps = [];
Cs_Xp_maps = [];

% read NPP data and convert into global map (unit: gC m-2 yr-1)
file_NPP = 'G:\My research\case2\working\cable-C\Matrix_out\';
% read ecosystem residence time data (unit: days)
file_tuaE = 'G:\My research\case2\working\cable-C\Transit_Matrix\';
% read C stroage capacity (unit: gC m-2)
file_Xc = 'G:\My research\case2\working\cable-C\Transit_Matrix\';
% read Xp data (unit: gC m-2)
file_Xp = 'G:\My research\case2\working\cable-C\Transit_Matrix\'
ndays = 365;


for i = 1901:2013
    i
    
    C_NPPinf = xlsread([file_NPP, 'Carbon_1_',num2str(i),'_run_CN_cycle.csv']);  
    C_Xcinf = xlsread([file_Xc,'C_capacity_',num2str(i),'.csv']);
    C_tuaE9 = xlsread([file_tuaE,'residet_',num2str(i),'.csv']);  % read residence time for each pools
    C_Xpinf = xlsread([file_Xp,'C_potencial_',num2str(i),'.csv']);
     
    Cp_NPP = C_NPPinf(:,1);             % C input into Plant C pools,         unit: gC m-2 yr-1
    Cp_Xc = nansum(C_Xcinf(:,1:3),2);   % plant C storage capacity,           unit: gC m-2
    Cp_tuaE = Cp_Xc./Cp_NPP;            % C residence time for plant C pools, unit: yr
    Cp_Xp = nansum(C_Xpinf(:,1:3),2);   % plant C storage potential,                unit: gC m-2
    
    Cs_Xc = nansum(C_Xcinf(:,4:9),2);   % soil C storage capacity, including litter and SOM, unit: gC m-2
    Cs_tuaE = nansum(C_tuaE9(:,4:9),2); % soil C residence time,                             unit: day
    Cs_tuaE = Cs_tuaE./ndays;           % covert unit into year,                             unit: year
    Cs_inC = Cs_Xc./Cs_tuaE;            % calculate C inputs into soil,                      unit: gC m-2 yr-1
    Cs_Xp = nansum(C_Xpinf(:,4:9),2);   % soil C storage potential,                          unit: gC m-2
    
    % convert into global map based on mask
    mask_Cp_NPP = mask; index_NPP_p = 1;
    mask_Cp_tuaE = mask; index_tuaE_p = 1;
    mask_Cp_Xc = mask; index_Xc_p = 1;
    mask_Cp_Xp = mask; index_Xp_p = 1;
    
    mask_Cs_inC = mask; index_NPP_s = 1;
    mask_Cs_tuaE = mask; index_tuaE_s = 1;
    mask_Cs_Xc = mask; index_Xc_s = 1;
    mask_Cs_Xp = mask; index_Xp_s = 1;
    
    for rowID=1:360;
        for colID=1:720
            % Plant
            if isnan(mask_Cp_NPP(rowID,colID)) == 0
               mask_Cp_NPP(rowID,colID) = Cp_NPP(index_NPP_p,1);
               index_NPP_p = index_NPP_p+1;
               index_NPP_p
            end
            
            if isnan(mask_Cp_tuaE(rowID,colID)) == 0
               mask_Cp_tuaE(rowID,colID) = Cp_tuaE(index_tuaE_p,1);
               index_tuaE_p = index_tuaE_p + 1;
            end
            
            if isnan(mask_Cp_Xc(rowID,colID)) == 0
               mask_Cp_Xc(rowID,colID) = Cp_Xc(index_Xc_p,1);
               index_Xc_p = index_Xc_p + 1;
            end
            
            if isnan(mask_Cp_Xp(rowID,colID)) == 0
               mask_Cp_Xp(rowID,colID) = Cp_Xp(index_Xp_p,1);
               index_Xp_p = index_Xp_p + 1;
            end
            
            % Soil
            if isnan(mask_Cs_inC(rowID,colID)) == 0
               mask_Cs_inC(rowID,colID) = Cs_inC(index_NPP_s,1);
               index_NPP_s = index_NPP_s+1;
               index_NPP_s
            end
            
            if isnan(mask_Cs_tuaE(rowID,colID)) == 0
               mask_Cs_tuaE(rowID,colID) = Cs_tuaE(index_tuaE_s,1);
               index_tuaE_s = index_tuaE_s + 1;
            end
            
            if isnan(mask_Cs_Xc(rowID,colID)) == 0
               mask_Cs_Xc(rowID,colID) = Cs_Xc(index_Xc_s,1);
               index_Xc_s = index_Xc_s + 1;
            end
            
            if isnan(mask_Cs_Xp(rowID,colID)) == 0
               mask_Cs_Xp(rowID,colID) = Cs_Xp(index_Xp_s,1);
               index_Xp_s = index_Xp_s + 1;
            end

        end
    end
    
    Cp_NPP_maps(:,:,i-1900) = mask_Cp_NPP;
    Cp_tuaE_maps(:,:,i-1900) = mask_Cp_tuaE;
    Cp_Xc_maps(:,:,i-1900) = mask_Cp_Xc;
    Cp_Xp_maps(:,:,i-1900) = mask_Cp_Xp;
    
    Cs_inC_maps(:,:,i-1900) = mask_Cs_inC;
    Cs_tuaE_maps(:,:,i-1900) = mask_Cs_tuaE;
    Cs_Xc_maps(:,:,i-1900) = mask_Cs_Xc;
    Cs_Xp_maps(:,:,i-1900) = mask_Cs_Xp;
end

clearvars -except mask cellarea ...
                  Cp_NPP_maps Cp_tuaE_maps Cp_Xc_maps Cp_Xp_maps Cs_inC_maps Cs_tuaE_maps Cs_Xc_maps Cs_Xp_maps
              
%% CN
% plant and soil
CNp_NPP_maps = [];
CNp_tuaE_maps = [];
CNp_Xc_maps = [];
CNp_Xp_maps = [];

CNs_inC_maps = [];
CNs_tuaE_maps = [];
CNs_Xc_maps = [];
CNs_Xp_maps = [];

% read NPP data and convert into global map (unit: gC m-2 yr-1)
file_NPP = 'G:\My research\case2\working\cable-CN\Matrix_out\';
% read ecosystem residence time data (unit: days)
file_tuaE = 'G:\My research\case2\working\cable-CN\Transit_Matrix\';
% read C stroage capacity (unit: gC m-2)
file_Xc = 'G:\My research\case2\working\cable-CN\Transit_Matrix\';
% read Xp data (unit: gC m-2)
file_Xp = 'G:\My research\case2\working\cable-CN\Transit_Matrix\';
ndays = 365;

for i = 1901:2013
    i
    
    CN_NPPinf = xlsread([file_NPP, 'Carbon_1_',num2str(i),'_run_CN_cycle.csv']);  
    CN_Xcinf = xlsread([file_Xc,'C_capacity_',num2str(i),'.csv']);
    CN_tuaE9 = xlsread([file_tuaE,'residet_',num2str(i),'.csv']);  % read residence time for each pools
    CN_Xpinf = xlsread([file_Xp,'C_potencial_',num2str(i),'.csv']);
     
    CNp_NPP = CN_NPPinf(:,1);             % C input into Plant C pools,         unit: gC m-2 yr-1
    CNp_Xc = nansum(CN_Xcinf(:,1:3),2);   % plant C storage capacity,           unit: gC m-2
    CNp_tuaE = CNp_Xc./CNp_NPP;           % C residence time for plant C pools, unit: yr
    CNp_Xp = nansum(CN_Xpinf(:,1:3),2);   % plant C storage potential,                unit: gC m-2
    
    CNs_Xc = nansum(CN_Xcinf(:,4:9),2);   % soil C storage capacity, including litter and SOM, unit: gC m-2
    CNs_tuaE = nansum(CN_tuaE9(:,4:9),2); % soil C residence time,                             unit: day
    CNs_tuaE = CNs_tuaE./ndays;           % covert unit into year,                             unit: year
    CNs_inC = CNs_Xc./CNs_tuaE;           % calculate C inputs into soil,                      unit: gC m-2 yr-1
    CNs_Xp = nansum(CN_Xpinf(:,4:9),2);   % soil C storage potential,                          unit: gC m-2
    
    % convert into global map based on mask
    mask_CNp_NPP = mask; index_NPP_p = 1;
    mask_CNp_tuaE = mask; index_tuaE_p = 1;
    mask_CNp_Xc = mask; index_Xc_p = 1;
    mask_CNp_Xp = mask; index_Xp_p = 1;
    
    mask_CNs_inC = mask; index_NPP_s = 1;
    mask_CNs_tuaE = mask; index_tuaE_s = 1;
    mask_CNs_Xc = mask; index_Xc_s = 1;
    mask_CNs_Xp = mask; index_Xp_s = 1;
    
    for rowID=1:360;
        for colID=1:720
            % Plant
            if isnan(mask_CNp_NPP(rowID,colID)) == 0
               mask_CNp_NPP(rowID,colID) = CNp_NPP(index_NPP_p,1);
               index_NPP_p = index_NPP_p+1;
               index_NPP_p
            end
            
            if isnan(mask_CNp_tuaE(rowID,colID)) == 0
               mask_CNp_tuaE(rowID,colID) = CNp_tuaE(index_tuaE_p,1);
               index_tuaE_p = index_tuaE_p + 1;
            end
            
            if isnan(mask_CNp_Xc(rowID,colID)) == 0
               mask_CNp_Xc(rowID,colID) = CNp_Xc(index_Xc_p,1);
               index_Xc_p = index_Xc_p + 1;
            end
            
            if isnan(mask_CNp_Xp(rowID,colID)) == 0
               mask_CNp_Xp(rowID,colID) = CNp_Xp(index_Xp_p,1);
               index_Xp_p = index_Xp_p + 1;
            end
            
            % Soil
            if isnan(mask_CNs_inC(rowID,colID)) == 0
               mask_CNs_inC(rowID,colID) = CNs_inC(index_NPP_s,1);
               index_NPP_s = index_NPP_s+1;
               
            end
            
            if isnan(mask_CNs_tuaE(rowID,colID)) == 0
               mask_CNs_tuaE(rowID,colID) = CNs_tuaE(index_tuaE_s,1);
               index_tuaE_s = index_tuaE_s + 1;
            end
            
            if isnan(mask_CNs_Xc(rowID,colID)) == 0
               mask_CNs_Xc(rowID,colID) = CNs_Xc(index_Xc_s,1);
               index_Xc_s = index_Xc_s + 1;
            end
            
            if isnan(mask_CNs_Xp(rowID,colID)) == 0
               mask_CNs_Xp(rowID,colID) = CNs_Xp(index_Xp_s,1);
               index_Xp_s = index_Xp_s + 1;
            end

        end
    end
    
    CNp_NPP_maps(:,:,i-1900) = mask_CNp_NPP;
    CNp_tuaE_maps(:,:,i-1900) = mask_CNp_tuaE;
    CNp_Xc_maps(:,:,i-1900) = mask_CNp_Xc;
    CNp_Xp_maps(:,:,i-1900) = mask_CNp_Xp;
    
    CNs_inC_maps(:,:,i-1900) = mask_CNs_inC;
    CNs_tuaE_maps(:,:,i-1900) = mask_CNs_tuaE;
    CNs_Xc_maps(:,:,i-1900) = mask_CNs_Xc;
    CNs_Xp_maps(:,:,i-1900) = mask_CNs_Xp;
    
end

clearvars -except mask cellarea ...
                  Cp_NPP_maps Cp_tuaE_maps Cp_Xc_maps Cp_Xp_maps Cs_inC_maps Cs_tuaE_maps Cs_Xc_maps Cs_Xp_maps ...
                  CNp_NPP_maps CNp_tuaE_maps CNp_Xc_maps CNp_Xp_maps CNs_inC_maps CNs_tuaE_maps CNs_Xc_maps CNs_Xp_maps
              
%% CNP
% plant and soil
CNPp_NPP_maps = [];
CNPp_tuaE_maps = [];
CNPp_Xc_maps = [];
CNPp_Xp_maps = [];

CNPs_inC_maps = [];
CNPs_tuaE_maps = [];
CNPs_Xc_maps = [];
CNPs_Xp_maps = [];

% read NPP data and convert into global map (unit: gC m-2 yr-1)
file_NPP = 'G:\My research\case2\working\cable-CNP\Matrix_out\';
% read ecosystem residence time data (unit: days)
file_tuaE = 'G:\My research\case2\working\cable-CNP\Transit_Matrix\';
% read C stroage capacity (unit: gC m-2)
file_Xc = 'G:\My research\case2\working\cable-CNP\Transit_Matrix\';
% read Xp data (unit: gC m-2)
file_Xp = 'G:\My research\case2\working\cable-CNP\Transit_Matrix\';

ndays = 365;

for i = 1901:2013
    i
    
    CNP_NPPinf = xlsread([file_NPP, 'Carbon_1_',num2str(i),'_run_CN_cycle.csv']);  
    CNP_Xcinf = xlsread([file_Xc,'C_capacity_',num2str(i),'.csv']);
    CNP_tuaE9 = xlsread([file_tuaE,'residet_',num2str(i),'.csv']);  % read residence time for each pools
    CNP_Xpinf = xlsread([file_Xp,'C_potencial_',num2str(i),'.csv']);
     
    CNPp_NPP = CNP_NPPinf(:,1);             % C input into Plant C pools,         unit: gC m-2 yr-1
    CNPp_Xc = nansum(CNP_Xcinf(:,1:3),2);   % plant C storage capacity,           unit: gC m-2
    CNPp_tuaE = CNPp_Xc./CNPp_NPP;          % C residence time for plant C pools, unit: yr
    CNPp_Xp = nansum(CNP_Xpinf(:,1:3),2);   % plant C storage potential,                unit: gC m-2
    
    CNPs_Xc = nansum(CNP_Xcinf(:,4:9),2);   % soil C storage capacity, including litter and SOM, unit: gC m-2
    CNPs_tuaE = nansum(CNP_tuaE9(:,4:9),2); % soil C residence time,                             unit: day
    CNPs_tuaE = CNPs_tuaE./ndays;           % covert unit into year,                             unit: year
    CNPs_inC = CNPs_Xc./CNPs_tuaE;          % calculate C inputs into soil,                      unit: gC m-2 yr-1
    CNPs_Xp = nansum(CNP_Xpinf(:,4:9),2);   % soil C storage potential,                          unit: gC m-2
    
    % convert into global map based on mask
    mask_CNPp_NPP = mask; index_NPP_p = 1;
    mask_CNPp_tuaE = mask; index_tuaE_p = 1;
    mask_CNPp_Xc = mask; index_Xc_p = 1;
    mask_CNPp_Xp = mask; index_Xp_p = 1;
    
    mask_CNPs_inC = mask; index_NPP_s = 1;
    mask_CNPs_tuaE = mask; index_tuaE_s = 1;
    mask_CNPs_Xc = mask; index_Xc_s = 1;
    mask_CNPs_Xp = mask; index_Xp_s = 1;
    
    for rowID=1:360;
        for colID=1:720
            % Plant
            if isnan(mask_CNPp_NPP(rowID,colID)) == 0
               mask_CNPp_NPP(rowID,colID) = CNPp_NPP(index_NPP_p,1);
               index_NPP_p = index_NPP_p+1;
               index_NPP_p
            end
            
            if isnan(mask_CNPp_tuaE(rowID,colID)) == 0
               mask_CNPp_tuaE(rowID,colID) = CNPp_tuaE(index_tuaE_p,1);
               index_tuaE_p = index_tuaE_p + 1;
            end
            
            if isnan(mask_CNPp_Xc(rowID,colID)) == 0
               mask_CNPp_Xc(rowID,colID) = CNPp_Xc(index_Xc_p,1);
               index_Xc_p = index_Xc_p + 1;
            end
            
            if isnan(mask_CNPp_Xp(rowID,colID)) == 0
               mask_CNPp_Xp(rowID,colID) = CNPp_Xp(index_Xp_p,1);
               index_Xp_p = index_Xp_p + 1;
            end
            
            % Soil
            if isnan(mask_CNPs_inC(rowID,colID)) == 0
               mask_CNPs_inC(rowID,colID) = CNPs_inC(index_NPP_s,1);
               index_NPP_s = index_NPP_s+1;
               
            end
            
            if isnan(mask_CNPs_tuaE(rowID,colID)) == 0
               mask_CNPs_tuaE(rowID,colID) = CNPs_tuaE(index_tuaE_s,1);
               index_tuaE_s = index_tuaE_s + 1;
            end
            
            if isnan(mask_CNPs_Xc(rowID,colID)) == 0
               mask_CNPs_Xc(rowID,colID) = CNPs_Xc(index_Xc_s,1);
               index_Xc_s = index_Xc_s + 1;
            end
            
            if isnan(mask_CNPs_Xp(rowID,colID)) == 0
               mask_CNPs_Xp(rowID,colID) = CNPs_Xp(index_Xp_s,1);
               index_Xp_s = index_Xp_s + 1;
            end

        end
    end
    
    CNPp_NPP_maps(:,:,i-1900) = mask_CNPp_NPP;
    CNPp_tuaE_maps(:,:,i-1900) = mask_CNPp_tuaE;
    CNPp_Xc_maps(:,:,i-1900) = mask_CNPp_Xc;
    CNPp_Xp_maps(:,:,i-1900) = mask_CNPp_Xp;
    
    CNPs_inC_maps(:,:,i-1900) = mask_CNPs_inC;
    CNPs_tuaE_maps(:,:,i-1900) = mask_CNPs_tuaE;
    CNPs_Xc_maps(:,:,i-1900) = mask_CNPs_Xc;
    CNPs_Xp_maps(:,:,i-1900) = mask_CNPs_Xp;
end

clearvars -except mask cellarea ...
                  Cp_NPP_maps Cp_tuaE_maps Cp_Xc_maps Cp_Xp_maps Cs_inC_maps Cs_tuaE_maps Cs_Xc_maps Cs_Xp_maps ...
                  CNp_NPP_maps CNp_tuaE_maps CNp_Xc_maps CNp_Xp_maps CNs_inC_maps CNs_tuaE_maps CNs_Xc_maps CNs_Xp_maps ...
                  CNPp_NPP_maps CNPp_tuaE_maps CNPp_Xc_maps CNPp_Xp_maps CNPs_inC_maps CNPs_tuaE_maps CNPs_Xc_maps CNPs_Xp_maps 

cd('G:\My research\case2\JGR\Working\step4_PlantSoil')   
save('CNPps_NPPtuaE_maps.mat')

clear
clc
load('G:\My research\case2\JGR\Working\step4_PlantSoil\CNPps_NPPtuaE_maps.mat')   
%% Figure
% calculate iniitial values of C inputs into plant and soil
% also Net changes

% C-only
% plant 
Cp_NPP_10bg = Cp_NPP_maps(:,:,1:10);
Cp_NPP_10ed = Cp_NPP_maps(:,:,104:113);
Cp_tuaE_10bg = Cp_tuaE_maps(:,:,1:10);
Cp_tuaE_10ed = Cp_tuaE_maps(:,:,104:113);
Cp_Xp_10bg = Cp_Xp_maps(:,:,1:10);
Cp_Xp_10ed = Cp_Xp_maps(:,:,104:113);

Cp_NPP_10bg_ag = nanmean(Cp_NPP_10bg,3);
Cp_NPP_10ed_ag = nanmean(Cp_NPP_10ed,3);
Cp_tuaE_10bg_ag = nanmean(Cp_tuaE_10bg,3);
Cp_tuaE_10ed_ag = nanmean(Cp_tuaE_10ed,3);
Cp_Xp_10bg_ag = nanmean(Cp_Xp_10bg,3);
Cp_Xp_10ed_ag = nanmean(Cp_Xp_10ed,3);

Net_Cp_NPP_map = Cp_NPP_10ed_ag - Cp_NPP_10bg_ag;         % Net change in plant C inputs,                    unit: gC m-2 yr-1
Net_Cp_tuaE_map = Cp_tuaE_10ed_ag - Cp_tuaE_10bg_ag;      % Net change in C residence time of plant C pool, unit: year
Net_Cp_Xp_map = Cp_Xp_10ed_ag - Cp_Xp_10bg_ag;            % Net change in plant Xp,                          unit: gC m-2 yr-1

% soil
Cs_inC_10bg = Cs_inC_maps(:,:,1:10);
Cs_inC_10ed = Cs_inC_maps(:,:,104:113);
Cs_tuaE_10bg = Cs_tuaE_maps(:,:,1:10);
Cs_tuaE_10ed = Cs_tuaE_maps(:,:,104:113);
Cs_Xp_10bg = Cs_Xp_maps(:,:,1:10);
Cs_Xp_10ed = Cs_Xp_maps(:,:,104:113);

Cs_inC_10bg_ag = nanmean(Cs_inC_10bg,3);
Cs_inC_10ed_ag = nanmean(Cs_inC_10ed,3);
Cs_tuaE_10bg_ag = nanmean(Cs_tuaE_10bg,3);
Cs_tuaE_10ed_ag = nanmean(Cs_tuaE_10ed,3);
Cs_Xp_10bg_ag = nanmean(Cs_Xp_10bg,3);
Cs_Xp_10ed_ag = nanmean(Cs_Xp_10ed,3);

Net_Cs_inC_map = Cs_inC_10ed_ag - Cs_inC_10bg_ag;         % Net change in soil C inputs,                    unit: gC m-2 yr-1 
Net_Cs_tuaE_map = Cs_tuaE_10ed_ag - Cs_tuaE_10bg_ag;      % Net change in C residence time of soil C pool, unit: year
Net_Cs_Xp_map = Cs_Xp_10ed_ag - Cs_Xp_10bg_ag;            % Net change in soil Xp,                          unit: gC m-2 yr-1

% CN
% plant 
CNp_NPP_10bg = CNp_NPP_maps(:,:,1:10);
CNp_NPP_10ed = CNp_NPP_maps(:,:,104:113);
CNp_tuaE_10bg = CNp_tuaE_maps(:,:,1:10);
CNp_tuaE_10ed = CNp_tuaE_maps(:,:,104:113);
CNp_Xp_10bg = CNp_Xp_maps(:,:,1:10);
CNp_Xp_10ed = CNp_Xp_maps(:,:,104:113);

CNp_NPP_10bg_ag = nanmean(CNp_NPP_10bg,3);
CNp_NPP_10ed_ag = nanmean(CNp_NPP_10ed,3);
CNp_tuaE_10bg_ag = nanmean(CNp_tuaE_10bg,3);
CNp_tuaE_10ed_ag = nanmean(CNp_tuaE_10ed,3);
CNp_Xp_10bg_ag = nanmean(CNp_Xp_10bg,3);
CNp_Xp_10ed_ag = nanmean(CNp_Xp_10ed,3);

Net_CNp_NPP_map = CNp_NPP_10ed_ag - CNp_NPP_10bg_ag;         % Net change in plant C inputs,                    unit: gC m-2 yr-1
Net_CNp_tuaE_map = CNp_tuaE_10ed_ag - CNp_tuaE_10bg_ag;      % Net change in C residence time of plant C pool, unit: year
Net_CNp_Xp_map = CNp_Xp_10ed_ag - CNp_Xp_10bg_ag;            % Net change in plant Xp,                          unit: gC m-2 yr-1

% soil
CNs_inC_10bg = CNs_inC_maps(:,:,1:10);
CNs_inC_10ed = CNs_inC_maps(:,:,104:113);
CNs_tuaE_10bg = CNs_tuaE_maps(:,:,1:10);
CNs_tuaE_10ed = CNs_tuaE_maps(:,:,104:113);
CNs_Xp_10bg = CNs_Xp_maps(:,:,1:10);
CNs_Xp_10ed = CNs_Xp_maps(:,:,104:113);

CNs_inC_10bg_ag = nanmean(CNs_inC_10bg,3);
CNs_inC_10ed_ag = nanmean(CNs_inC_10ed,3);
CNs_tuaE_10bg_ag = nanmean(CNs_tuaE_10bg,3);
CNs_tuaE_10ed_ag = nanmean(CNs_tuaE_10ed,3);
CNs_Xp_10bg_ag = nanmean(CNs_Xp_10bg,3);
CNs_Xp_10ed_ag = nanmean(CNs_Xp_10ed,3);

Net_CNs_inC_map = CNs_inC_10ed_ag - CNs_inC_10bg_ag;         % Net change in soil C inputs,                    unit: gC m-2 yr-1 
Net_CNs_tuaE_map = CNs_tuaE_10ed_ag - CNs_tuaE_10bg_ag;      % Net change in C residence time of soil C pool, unit: year
Net_CNs_Xp_map = CNs_Xp_10ed_ag - CNs_Xp_10bg_ag;            % Net change in soil Xp,                          unit: gC m-2 yr-1

% CNP
% plant 
CNPp_NPP_10bg = CNPp_NPP_maps(:,:,1:10);
CNPp_NPP_10ed = CNPp_NPP_maps(:,:,104:113);
CNPp_tuaE_10bg = CNPp_tuaE_maps(:,:,1:10);
CNPp_tuaE_10ed = CNPp_tuaE_maps(:,:,104:113);
CNPp_Xp_10bg = CNPp_Xp_maps(:,:,1:10);
CNPp_Xp_10ed = CNPp_Xp_maps(:,:,104:113);

CNPp_NPP_10bg_ag = nanmean(CNPp_NPP_10bg,3);
CNPp_NPP_10ed_ag = nanmean(CNPp_NPP_10ed,3);
CNPp_tuaE_10bg_ag = nanmean(CNPp_tuaE_10bg,3);
CNPp_tuaE_10ed_ag = nanmean(CNPp_tuaE_10ed,3);
CNPp_Xp_10bg_ag = nanmean(CNPp_Xp_10bg,3);
CNPp_Xp_10ed_ag = nanmean(CNPp_Xp_10ed,3);

Net_CNPp_NPP_map = CNPp_NPP_10ed_ag - CNPp_NPP_10bg_ag;         % Net change in plant C inputs,                    unit: gC m-2 yr-1
Net_CNPp_tuaE_map = CNPp_tuaE_10ed_ag - CNPp_tuaE_10bg_ag;      % Net change in C residence time of plant C pool, unit: year
Net_CNPp_Xp_map = CNPp_Xp_10ed_ag - CNPp_Xp_10bg_ag;            % Net change in plant Xp,                          unit: gC m-2 yr-1

% soil
CNPs_inC_10bg = CNPs_inC_maps(:,:,1:10);
CNPs_inC_10ed = CNPs_inC_maps(:,:,104:113);
CNPs_tuaE_10bg = CNPs_tuaE_maps(:,:,1:10);
CNPs_tuaE_10ed = CNPs_tuaE_maps(:,:,104:113);
CNPs_Xp_10bg = CNPs_Xp_maps(:,:,1:10);
CNPs_Xp_10ed = CNPs_Xp_maps(:,:,104:113);

CNPs_inC_10bg_ag = nanmean(CNPs_inC_10bg,3);
CNPs_inC_10ed_ag = nanmean(CNPs_inC_10ed,3);
CNPs_tuaE_10bg_ag = nanmean(CNPs_tuaE_10bg,3);
CNPs_tuaE_10ed_ag = nanmean(CNPs_tuaE_10ed,3);
CNPs_Xp_10bg_ag = nanmean(CNPs_Xp_10bg,3);
CNPs_Xp_10ed_ag = nanmean(CNPs_Xp_10ed,3);

Net_CNPs_inC_map = CNPs_inC_10ed_ag - CNPs_inC_10bg_ag;         % Net change in soil C inputs,                    unit: gC m-2 yr-1 
Net_CNPs_tuaE_map = CNPs_tuaE_10ed_ag - CNPs_tuaE_10bg_ag;      % Net change in C residence time of soil C pool, unit: year
Net_CNPs_Xp_map = CNPs_Xp_10ed_ag - CNPs_Xp_10bg_ag;            % Net change in soil Xp,                          unit: gC m-2 yr-1

% Figure for initial valures of C input and residence time (C-only, CN and CNP)
% preparing for data

Data_inCtuaE0_ps(:,:,1) = Cp_NPP_10bg_ag;
Data_inCtuaE0_ps(:,:,2) = Cs_inC_10bg_ag;
Data_inCtuaE0_ps(:,:,3) = Cp_tuaE_10bg_ag;
Data_inCtuaE0_ps(:,:,4) = Cs_tuaE_10bg_ag;

Data_inCtuaE0_ps(:,:,5) = CNp_NPP_10bg_ag;
Data_inCtuaE0_ps(:,:,6) = CNs_inC_10bg_ag;
Data_inCtuaE0_ps(:,:,7) = CNp_tuaE_10bg_ag;
Data_inCtuaE0_ps(:,:,8) = CNs_tuaE_10bg_ag;

Data_inCtuaE0_ps(:,:,9) = CNPp_NPP_10bg_ag;
Data_inCtuaE0_ps(:,:,10) = CNPs_inC_10bg_ag;
Data_inCtuaE0_ps(:,:,11) = CNPp_tuaE_10bg_ag;
Data_inCtuaE0_ps(:,:,12) = CNPs_tuaE_10bg_ag;

Data_inCtuaE0_ps(301:360,:,:) = [];

%ax = gca;
%mymap_inC0_sp = colormap(ax)
%save('mymap_inC0_sp','mymap_inC0_sp')
cd('G:\My research\case2\JGR\Working\step4_PlantSoil\inCtuaE0_soilplant')
% global maps for initial values of C inputs and residence time (plant and soil) 
%load('F:\My research\case2\JGR\Working\step3_reNPPtuaE_fct\mymap_CNP_re.mat')
load('G:\My research\case2\JGR\Working\step4_PlantSoil\inCtuaE0_soilplant\mymap_inC0_sp.mat')
load('G:\My research\case2\JGR\Working\step2_NPP0_tuaE0\mymap_CNP_tuaE2.mat')

figure
set(gcf,'position',[100 100 950 450])
maps_CNP = tight_subplot(3,4,[-0.08 -0.08],[0.06 0.02],[0.01 0.01])

Labels = {'(a)','(b)','(c)','(d)',...
          '(e)','(f)','(g)','(h)',...
          '(i)','(j)','(k)','(l)',...
          '(m)','(n)','(o)','(p)',...
          '(q)','(r)','(s)','(t)'};

% C inputs into plant C pool (original unit: gC m-2 yr-1)
for i= 1:4:9
    i
    
    map_i = flipud(Data_inCtuaE0_ps(:,:,i))./1000;   % Unit: KgC m-2 yr-1
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_inC0_sp)
    caxis([0 2])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end      

% C inputs into soil C pool (original unit: gC m-2 yr-1)
for i= 2:4:10
    i
    
    map_i = flipud(Data_inCtuaE0_ps(:,:,i))./1000;   % Unit: KgC m-2 yr-1
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_inC0_sp)
    caxis([0 2])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

% C residence time in plant pool (unit: year)
for i= 3:4:11
    i
    
    map_i = flipud(Data_inCtuaE0_ps(:,:,i));
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_CNP_tuaE2)
    caxis([0 24])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

% C residence time in soil pool (unit: year)
for i= 4:4:12
    i
    
    map_i = flipud(Data_inCtuaE0_ps(:,:,i));
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_CNP_tuaE2)
    caxis([0 600])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end


%% Figure
% Nutrient limitation on initial Cin and residence time (plant and soil)

% plant
Nlimt_NPP0_p = CNp_NPP_10bg_ag - Cp_NPP_10bg_ag;
Nlimt_tuaE0_p = CNp_tuaE_10bg_ag - Cp_tuaE_10bg_ag;

Plimt_NPP0_p = CNPp_NPP_10bg_ag - CNp_NPP_10bg_ag;
Plimt_tuaE0_p = CNPp_tuaE_10bg_ag - CNp_tuaE_10bg_ag;

% soil
Nlimt_NPP0_s = CNs_inC_10bg_ag - Cs_inC_10bg_ag;
Nlimt_tuaE0_s = CNs_tuaE_10bg_ag - Cs_tuaE_10bg_ag;

Plimt_inC0_s = CNPs_inC_10bg_ag - CNs_inC_10bg_ag;
Plimt_tuaE0_s = CNPs_tuaE_10bg_ag - CNs_tuaE_10bg_ag;

% Plant
Data_NPlit_p(:,:,1) = Nlimt_NPP0_p;
Data_NPlit_p(:,:,2) = Nlimt_tuaE0_p;
Data_NPlit_p(:,:,3) = Plimt_NPP0_p;
Data_NPlit_p(:,:,4) = Plimt_tuaE0_p;

Data_NPlit_p(301:360,:,:) = [];

% Soil
Data_NPlit_s(:,:,1) = Nlimt_NPP0_s;
Data_NPlit_s(:,:,2) = Nlimt_tuaE0_s;
Data_NPlit_s(:,:,3) = Plimt_inC0_s;
Data_NPlit_s(:,:,4) = Plimt_tuaE0_s;

Data_NPlit_s(301:360,:,:) = [];

% Figure for plant
figure
set(gcf,'position',[100 100 650 420])
maps_CNP = tight_subplot(2,2,[-0.14 -0.09],[0.08 0.001],[0.02 0.02])
load 'G:\My research\case2\JGR\Working\step2_NPP0_tuaE0\mymap_NP_limt0.mat'

Labels = {'(a)','(b)','(c)','(d)',...
          '(e)','(f)','(g)','(h)',...
          '(i)','(j)','(k)','(l)',...
          '(m)','(n)','(o)','(p)',...
          '(q)','(r)','(s)','(t)'};

for i= [1 3]
    i
    
    map_i = flipud(Data_NPlit_p(:,:,i))./1000;   % Unit: KgC m-2 yr-1
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP,mymap_NP_limt0)
    caxis([-1.5 1.5])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end
      
for i= [2 4]
    i
    
    map_i = flipud(Data_NPlit_p(:,:,i));   % Unit: years
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP,mymap_NP_limt0)
    caxis([-5 5])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end


% Figure for soil
figure
set(gcf,'position',[100 100 650 420])
maps_CNP = tight_subplot(2,2,[-0.14 -0.09],[0.08 0.001],[0.02 0.02])
load 'F:\My research\case2\JGR\Working\step2_NPP0_tuaE0\mymap_NP_limt0.mat'

Labels = {'(a)','(b)','(c)','(d)',...
          '(e)','(f)','(g)','(h)',...
          '(i)','(j)','(k)','(l)',...
          '(m)','(n)','(o)','(p)',...
          '(q)','(r)','(s)','(t)'};

for i= [1 3]
    i
    
    map_i = flipud(Data_NPlit_s(:,:,i))./1000;   % Unit: KgC m-2 yr-1
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP,mymap_NP_limt0)
    caxis([-1.5 1.5])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end
      
for i= [2 4]
    i
    
    map_i = flipud(Data_NPlit_s(:,:,i));   % Unit: years
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP,mymap_NP_limt0)
    caxis([-20 20])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

%% Figure
% re_NPP * tuaE0, re_tuaE * NPP0,  re_NPP * re_tuaE and Xp
% for plant and soil

% plant
% C-only
fct_reNPP_Cp = Net_Cp_NPP_map.* Cp_tuaE_10bg_ag;
fct_retuaE_Cp = Net_Cp_tuaE_map.* Cp_NPP_10bg_ag;
fct_reNPPtuaE_Cp = Net_Cp_NPP_map.* Net_Cp_tuaE_map;
fct_Xp_Cp = Net_Cp_Xp_map;

% CN
fct_reNPP_CNp = Net_CNp_NPP_map.* CNp_tuaE_10bg_ag;
fct_retuaE_CNp = Net_CNp_tuaE_map.* CNp_NPP_10bg_ag;
fct_reNPPtuaE_CNp = Net_CNp_NPP_map.* Net_CNp_tuaE_map;
fct_Xp_CNp = Net_CNp_Xp_map;

% CNP
fct_reNPP_CNPp = Net_CNPp_NPP_map.* CNPp_tuaE_10bg_ag;
fct_retuaE_CNPp = Net_CNPp_tuaE_map.* CNPp_NPP_10bg_ag;
fct_reNPPtuaE_CNPp = Net_CNPp_NPP_map.* Net_CNPp_tuaE_map;
fct_Xp_CNPp = Net_CNPp_Xp_map;

% preparing for dataset 
Data_reNPPtuaE_p(:,:,1) = fct_reNPP_Cp;
Data_reNPPtuaE_p(:,:,2) = fct_retuaE_Cp;
Data_reNPPtuaE_p(:,:,3) = fct_reNPPtuaE_Cp;
Data_reNPPtuaE_p(:,:,4) = fct_Xp_Cp;

Data_reNPPtuaE_p(:,:,5) = fct_reNPP_CNp;
Data_reNPPtuaE_p(:,:,6) = fct_retuaE_CNp;
Data_reNPPtuaE_p(:,:,7) = fct_reNPPtuaE_CNp;
Data_reNPPtuaE_p(:,:,8) = fct_Xp_CNp;

Data_reNPPtuaE_p(:,:,9) = fct_reNPP_CNPp;
Data_reNPPtuaE_p(:,:,10) = fct_retuaE_CNPp;
Data_reNPPtuaE_p(:,:,11) = fct_reNPPtuaE_CNPp;
Data_reNPPtuaE_p(:,:,12) = fct_Xp_CNPp;

Data_reNPPtuaE_p(301:360,:,:) = [];

% figure
cd('G:\My research\case2\JGR\Working\step4_PlantSoil\re_inCtua_soilplant')
load('G:\My research\case2\JGR\Working\step3_reNPPtuaE_fct\mymap_CNP_re.mat')

figure
set(gcf,'position',[100 100 950 450])
maps_CNP = tight_subplot(3,4,[-0.08 -0.08],[0.06 0.02],[0.01 0.01])

Labels = {'(a)','(b)','(c)','(d)',...
          '(e)','(f)','(g)','(h)',...
          '(i)','(j)','(k)','(l)',...
          '(m)','(n)','(o)','(p)',...
          '(q)','(r)','(s)','(t)'};
      
for i= 1:4:9
    i
    
    map_i = flipud(Data_reNPPtuaE_p(:,:,i))./1000;   % Unit: KgC m-2 
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_CNP_re)
    caxis([-10 10])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end


for i= 2:4:10
    i
    
    map_i = flipud(Data_reNPPtuaE_p(:,:,i))./1000;   % Unit: KgC m-2
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_CNP_re)
    caxis([-10, 10])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end


for i= 3:4:11
    i
    
    map_i = flipud(Data_reNPPtuaE_p(:,:,i))./1000;   % Unit: KgC m-2
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_CNP_re)
    caxis([-0.2 0.2])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

for i= 4:4:12
    i
    
    map_i = flipud(Data_reNPPtuaE_p(:,:,i))./1000;   % Unit: KgC m-2
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_CNP_re)
    caxis([-10 10])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end


% Soil
% C-only
fct_reNPP_Cs = Net_Cs_inC_map.* Cs_tuaE_10bg_ag;
fct_retuaE_Cs = Net_Cs_tuaE_map.* Cs_inC_10bg_ag;
fct_reNPPtuaE_Cs = Net_Cs_inC_map.* Net_Cs_tuaE_map;
fct_Xp_Cs = Net_Cs_Xp_map;

% CN
fct_reNPP_CNs = Net_CNs_inC_map.* CNs_tuaE_10bg_ag;
fct_retuaE_CNs = Net_CNs_tuaE_map.* CNs_inC_10bg_ag;
fct_reNPPtuaE_CNs = Net_CNs_inC_map.* Net_CNs_tuaE_map;
fct_Xp_CNs = Net_CNs_Xp_map;

% CNP
fct_reNPP_CNPs = Net_CNPs_inC_map.* CNPs_tuaE_10bg_ag;
fct_retuaE_CNPs = Net_CNPs_tuaE_map.* CNPs_inC_10bg_ag;
fct_reNPPtuaE_CNPs = Net_CNPs_inC_map.* Net_CNPs_tuaE_map;
fct_Xp_CNPs = Net_CNPs_Xp_map;

% preparing for dataset 
Data_reNPPtuaE_s(:,:,1) = fct_reNPP_Cs;
Data_reNPPtuaE_s(:,:,2) = fct_retuaE_Cs;
Data_reNPPtuaE_s(:,:,3) = fct_reNPPtuaE_Cs;
Data_reNPPtuaE_s(:,:,4) = fct_Xp_Cs;

Data_reNPPtuaE_s(:,:,5) = fct_reNPP_CNs;
Data_reNPPtuaE_s(:,:,6) = fct_retuaE_CNs;
Data_reNPPtuaE_s(:,:,7) = fct_reNPPtuaE_CNs;
Data_reNPPtuaE_s(:,:,8) = fct_Xp_CNs;

Data_reNPPtuaE_s(:,:,9) = fct_reNPP_CNPs;
Data_reNPPtuaE_s(:,:,10) = fct_retuaE_CNPs;
Data_reNPPtuaE_s(:,:,11) = fct_reNPPtuaE_CNPs;
Data_reNPPtuaE_s(:,:,12) = fct_Xp_CNPs;

Data_reNPPtuaE_s(301:360,:,:) = [];

% figure
cd('F:\My research\case2\JGR\Working\step4_PlantSoil\re_inCtua_soilplant')
load('F:\My research\case2\JGR\Working\step3_reNPPtuaE_fct\mymap_CNP_re.mat')

figure
set(gcf,'position',[100 100 950 450])
maps_CNP = tight_subplot(3,4,[-0.08 -0.08],[0.06 0.02],[0.01 0.01])

Labels = {'(a)','(b)','(c)','(d)',...
          '(e)','(f)','(g)','(h)',...
          '(i)','(j)','(k)','(l)',...
          '(m)','(n)','(o)','(p)',...
          '(q)','(r)','(s)','(t)'};
      
for i= 1:4:9
    i
    
    map_i = flipud(Data_reNPPtuaE_s(:,:,i))./1000;   % Unit: KgC m-2 
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_CNP_re)
    caxis([-10 10])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end


for i= 2:4:10
    i
    
    map_i = flipud(Data_reNPPtuaE_s(:,:,i))./1000;   % Unit: KgC m-2
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_CNP_re)
    caxis([-0.2, 0.2])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end


for i= 3:4:11
    i
    
    map_i = flipud(Data_reNPPtuaE_s(:,:,i))./1000;   % Unit: KgC m-2
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_CNP_re)
    caxis([-0.2 0.2])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

for i= 4:4:12
    i
    
    map_i = flipud(Data_reNPPtuaE_s(:,:,i))./1000;   % Unit: KgC m-2
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_CNP_re)
    caxis([-10 10])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

%% Figure for Net change in NPP and residence time under three different coupling schemes (plant and soil)

% plant
% preparing for data
DataMap_NPPtuaE_reP(:,:,1) = Net_Cp_NPP_map;
DataMap_NPPtuaE_reP(:,:,2) = Net_Cp_tuaE_map;

DataMap_NPPtuaE_reP(:,:,3) = Net_CNp_NPP_map;
DataMap_NPPtuaE_reP(:,:,4) = Net_CNp_tuaE_map;

DataMap_NPPtuaE_reP(:,:,5) = Net_CNPp_NPP_map;
DataMap_NPPtuaE_reP(:,:,6) = Net_CNPp_tuaE_map;

DataMap_NPPtuaE_reP(301:360,:,:) = [];

load 'F:\My research\case2\JGR\Working\step2_NPP0_tuaE0\mymap_NP_limt0.mat'
figure
set(gcf,'position',[100 100 650 600])
maps_CNP = tight_subplot(3,2,[-0.08 -0.09],[0.06 0.02],[0.01 0.01])

Labels = {'(a)','(b)','(c)','(d)',...
          '(e)','(f)','(g)','(h)',...
          '(i)','(j)','(k)','(l)',...
          '(m)','(n)','(o)','(p)',...
          '(q)','(r)','(s)','(t)'};
      
for i= 1:2:5
    i
    
    map_i = flipud(DataMap_NPPtuaE_reP(:,:,i))./1000;   % Unit: KgC m-2 yr-1
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_NP_limt0)
    caxis([-0.5 0.5])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

for i= 2:2:6
    i
    
    map_i = flipud(DataMap_NPPtuaE_reP(:,:,i));   % Unit: year
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP,mymap_NP_limt0)
    caxis([-3 3])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

% soil
% preparing for data
DataMap_NPPtuaE_reS(:,:,1) = Net_Cs_inC_map;
DataMap_NPPtuaE_reS(:,:,2) = Net_Cs_tuaE_map;

DataMap_NPPtuaE_reS(:,:,3) = Net_CNs_inC_map;
DataMap_NPPtuaE_reS(:,:,4) = Net_CNs_tuaE_map;

DataMap_NPPtuaE_reS(:,:,5) = Net_CNPs_inC_map;
DataMap_NPPtuaE_reS(:,:,6) = Net_CNPs_tuaE_map;

DataMap_NPPtuaE_reS(301:360,:,:) = [];

figure
set(gcf,'position',[100 100 650 600])
maps_CNP = tight_subplot(3,2,[-0.08 -0.09],[0.06 0.02],[0.01 0.01])

Labels = {'(a)','(b)','(c)','(d)',...
          '(e)','(f)','(g)','(h)',...
          '(i)','(j)','(k)','(l)',...
          '(m)','(n)','(o)','(p)',...
          '(q)','(r)','(s)','(t)'};
      
for i= 1:2:5
    i
    
    map_i = flipud(DataMap_NPPtuaE_reS(:,:,i))./1000;   % Unit: KgC m-2 yr-1
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP, mymap_NP_limt0)
    caxis([-0.5 0.5])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

for i= 2:2:6
    i
    
    map_i = flipud(DataMap_NPPtuaE_reS(:,:,i));   % Unit: year
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    figCNP = maps_CNP(i);
    axes(figCNP);
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(figCNP,mymap_NP_limt0)
    caxis([-3 3])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end


















