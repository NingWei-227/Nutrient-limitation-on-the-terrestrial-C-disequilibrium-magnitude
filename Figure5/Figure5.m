clear
clc

% codes for global maps of plant C sinka and Soil C sink
% Comparison of land C sink under C-only, C-N and C-N-P interations
% read mask and cellarea data used to convert model outputs into global map
[num,text,raw]=xlsread('F:\My research\case2\working\ALL\codes\yanmo0.5.csv');
for rowID=1:360
    for colID =1:720
        if raw{rowID,colID}=='NA'
            raw{rowID,colID}=-9999;
        end  
    end
end
mask = cell2mat(raw);
mask(mask==-9999) =  nan;
cellarea = xlsread('F:\My research\case2\working\ALL\codes\area0.5.csv');   % unit: km2 = 10^6 m2
clearvars -except mask cellarea

%% C-only 
cd('F:\My research\case2\working\cable-C\Transit_Matrix');

Cp_pol_maps = [];
Cs_pol_maps = [];

% read data of C storege and convert into global maps (unit: gC m-2 yr-1)
for i = 1901:2013
    i
    
    C_pol9 = xlsread(['X_xx_',num2str(i),'.csv']);         % C storage of 9 C pools
    
    % plant C pool and soil C pool
    Cp_polE = sum(C_pol9(:,1:3),2);    % leaf, root and wood                           
    Cs_polE = sum(C_pol9(:,4:9),2);    % 3 litter pools and 3 SOM pools
    
    
    % convert into global map based on mask
    mask_Cp_pol = mask; index_cap = 1;
    mask_Cs_pol = mask; index_pol = 1;
    
    
    for rowID=1:360;
        for colID=1:720
            
            if isnan(mask_Cp_pol(rowID,colID)) == 0
               mask_Cp_pol(rowID,colID) = Cp_polE(index_cap,1);
               index_cap = index_cap+1;
               index_cap
            end
            if isnan(mask_Cs_pol(rowID,colID)) == 0
               mask_Cs_pol(rowID,colID) = Cs_polE(index_pol,1);
               index_pol = index_pol + 1;
            end
            
        end
    end
    
    Cp_pol_maps(:,:,i-1900) = mask_Cp_pol;
    Cs_pol_maps(:,:,i-1900) = mask_Cs_pol;
      
end

% global estimation of C storage (plant and soil)
Cp_pol_yr = [];
Cs_pol_yr = [];

% plant C sink and soil C sink
Cp_pol_yr_sk = [];
Cs_pol_yr_sk =[];

for yr = 1:113
    yr 
    
    Cp_cap_gb=  Cp_pol_maps(:,:,yr).*cellarea .* 10^6./10^15;     % convert unit from gC m-2 into Pg C
    Cs_pot_gb = Cs_pol_maps(:,:,yr).*cellarea .* 10^6./10^15;
    
    Cp_pol_yr(yr)  = nansum(Cp_cap_gb(:));
    Cs_pol_yr(yr)  = nansum(Cs_pot_gb(:));
    
end

Cp_pol_yr_sk = nanmean(Cp_pol_yr(104:113)) - nanmean(Cp_pol_yr(1:10));
Cs_pol_yr_sk = nanmean(Cs_pol_yr(104:113)) - nanmean(Cs_pol_yr(1:10));

clearvars -except mask cellarea Cp_pol_maps Cs_pol_maps Cp_pol_yr Cs_pol_yr Cp_pol_yr_sk Cs_pol_yr_sk


%% CN
cd('F:\My research\case2\working\cable-CN\Transit_Matrix');

CNp_pol_maps = [];
CNs_pol_maps = [];

% read data of C storege and convert into global maps (unit: gC m-2 yr-1)
for i = 1901:2013
    i
    
    CN_pol9 = xlsread(['X_xx_',num2str(i),'.csv']);         % C storage of 9 C pools
    
    % plant C pool and soil C pool
    CNp_polE = sum(CN_pol9(:,1:3),2);    % leaf, root and wood                           
    CNs_polE = sum(CN_pol9(:,4:9),2);    % 3 litter pools and 3 SOM pools
    
    
    % convert into global map based on mask
    mask_CNp_pol = mask; index_cap = 1;
    mask_CNs_pol = mask; index_pol = 1;
    
    
    for rowID=1:360;
        for colID=1:720
            
            if isnan(mask_CNp_pol(rowID,colID)) == 0
               mask_CNp_pol(rowID,colID) = CNp_polE(index_cap,1);
               index_cap = index_cap+1;
               index_cap
            end
            if isnan(mask_CNs_pol(rowID,colID)) == 0
               mask_CNs_pol(rowID,colID) = CNs_polE(index_pol,1);
               index_pol = index_pol + 1;
            end
            
        end
    end
    
    CNp_pol_maps(:,:,i-1900) = mask_CNp_pol;
    CNs_pol_maps(:,:,i-1900) = mask_CNs_pol;
      
end

% global estimation of C storage (plant and soil)
CNp_pol_yr = [];
CNs_pol_yr = [];

% plant C sink and soil C sink
CNp_pol_yr_sk = [];
CNs_pol_yr_sk =[];

for yr = 1:113
    yr 
    
    CNp_cap_gb=  CNp_pol_maps(:,:,yr).*cellarea .* 10^6./10^15;     % convert unit from gC m-2 into Pg C
    CNs_pot_gb = CNs_pol_maps(:,:,yr).*cellarea .* 10^6./10^15;
    
    CNp_pol_yr(yr)  = nansum(CNp_cap_gb(:));
    CNs_pol_yr(yr)  = nansum(CNs_pot_gb(:));
    
end

CNp_pol_yr_sk = nanmean(CNp_pol_yr(104:113)) - nanmean(CNp_pol_yr(1:10));
CNs_pol_yr_sk = nanmean(CNs_pol_yr(104:113)) - nanmean(CNs_pol_yr(1:10));

clearvars -except mask cellarea Cp_pol_maps Cs_pol_maps Cp_pol_yr Cs_pol_yr Cp_pol_yr_sk Cs_pol_yr_sk ...
                  CNp_pol_maps CNs_pol_maps CNp_pol_yr CNs_pol_yr CNp_pol_yr_sk CNs_pol_yr_sk 
save('F:\My research\case2\JGR\1_Figures\Figure5\Plant_Soil_X_Xc_Xp.mat') 
%% CNP
cd('F:\My research\case2\working\cable-CNP\Transit_Matrix');

CNPp_pol_maps = [];
CNPs_pol_maps = [];

% read data of C storege and convert into global maps (unit: gC m-2 yr-1)
for i = 1901:2013
    i
    
    CNP_pol9 = xlsread(['X_xx_',num2str(i),'.csv']);         % C storage of 9 C pools
    
    % plant C pool and soil C pool
    CNPp_polE = sum(CNP_pol9(:,1:3),2);    % leaf, root and wood                           
    CNPs_polE = sum(CNP_pol9(:,4:9),2);    % 3 litter pools and 3 SOM pools
    
    
    % convert into global map based on mask
    mask_CNPp_pol = mask; index_cap = 1;
    mask_CNPs_pol = mask; index_pol = 1;
    
    
    for rowID=1:360;
        for colID=1:720
            
            if isnan(mask_CNPp_pol(rowID,colID)) == 0
               mask_CNPp_pol(rowID,colID) = CNPp_polE(index_cap,1);
               index_cap = index_cap+1;
               index_cap
            end
            if isnan(mask_CNPs_pol(rowID,colID)) == 0
               mask_CNPs_pol(rowID,colID) = CNPs_polE(index_pol,1);
               index_pol = index_pol + 1;
            end
            
        end
    end
    
    CNPp_pol_maps(:,:,i-1900) = mask_CNPp_pol;
    CNPs_pol_maps(:,:,i-1900) = mask_CNPs_pol;
      
end

% global estimation of C storage (plant and soil)
CNPp_pol_yr = [];
CNPs_pol_yr = [];

% plant C sink and soil C sink
CNPp_pol_yr_sk = [];
CNPs_pol_yr_sk =[];

for yr = 1:113
    yr 
    
    CNPp_cap_gb=  CNPp_pol_maps(:,:,yr).*cellarea .* 10^6./10^15;     % convert unit from gC m-2 into Pg C
    CNPs_pot_gb = CNPs_pol_maps(:,:,yr).*cellarea .* 10^6./10^15;
    
    CNPp_pol_yr(yr)  = nansum(CNPp_cap_gb(:));
    CNPs_pol_yr(yr)  = nansum(CNPs_pot_gb(:));
    
end

CNPp_pol_yr_sk = nanmean(CNPp_pol_yr(104:113)) - nanmean(CNPp_pol_yr(1:10));
CNPs_pol_yr_sk = nanmean(CNPs_pol_yr(104:113)) - nanmean(CNPs_pol_yr(1:10));

clearvars -except mask cellarea Cp_pol_maps Cs_pol_maps Cp_pol_yr Cs_pol_yr Cp_pol_yr_sk Cs_pol_yr_sk ...
                  CNp_pol_maps CNs_pol_maps CNp_pol_yr CNs_pol_yr CNp_pol_yr_sk CNs_pol_yr_sk ...
                  CNPp_pol_maps CNPs_pol_maps CNPp_pol_yr CNPs_pol_yr CNPp_pol_yr_sk CNPs_pol_yr_sk  
load('F:\My research\case2\JGR\Working\step4_PlantSoil\CNPps_NPPtuaE_maps.mat') 

clearvars -except mask cellarea Cp_pol_maps Cs_pol_maps Cp_pol_yr Cs_pol_yr Cp_pol_yr_sk Cs_pol_yr_sk ...
                  CNp_pol_maps CNs_pol_maps CNp_pol_yr CNs_pol_yr CNp_pol_yr_sk CNs_pol_yr_sk ...
                  CNPp_pol_maps CNPs_pol_maps CNPp_pol_yr CNPs_pol_yr CNPp_pol_yr_sk CNPs_pol_yr_sk ...
                  Cp_Xp_maps Cs_Xp_maps Cp_Xp_yr Cs_Xp_yr Cp_Xp_yr_sk Cs_Xp_yr_sk ...
                  CNp_Xp_maps CNs_Xp_maps CNp_Xp_yr CNs_Xp_yr CNp_Xp_yr_sk CNs_Xp_yr_sk ...
                  CNPp_Xp_maps CNPs_Xp_maps CNPp_Xp_yr CNPs_Xp_yr CNPp_Xp_yr_sk CNPs_Xp_yr_sk


save('F:\My research\case2\JGR\1_Figures\Figure5\Plant_Soil_X_Xc_Xp.mat')              
%% Figure5:
clear;clc
load('F:\My research\case2\JGR\1_Figures\Figure5\Plant_Soil_X_Xc_Xp.mat')
% the spatial distribution of the magnitude of disequilibrium for plant and soil
% C-only
Cp_Xp_sp_sk = nanmean(Cp_Xp_maps(:,:,104:113),3) - nanmean(Cp_Xp_maps(:,:,1:10),3);
Cs_Xp_sp_sk = nanmean(Cs_Xp_maps(:,:,104:113),3) - nanmean(Cs_Xp_maps(:,:,1:10),3);
% CN
CNp_Xp_sp_sk = nanmean(CNp_Xp_maps(:,:,104:113),3) - nanmean(CNp_Xp_maps(:,:,1:10),3);
CNs_Xp_sp_sk = nanmean(CNs_Xp_maps(:,:,104:113),3) - nanmean(CNs_Xp_maps(:,:,1:10),3);
% CNP
CNPp_Xp_sp_sk = nanmean(CNPp_Xp_maps(:,:,104:113),3) - nanmean(CNPp_Xp_maps(:,:,1:10),3);
CNPs_Xp_sp_sk = nanmean(CNPs_Xp_maps(:,:,104:113),3) - nanmean(CNPs_Xp_maps(:,:,1:10),3);

Data_Xp_ps_map(:,:,1) = Cp_Xp_sp_sk;
Data_Xp_ps_map(:,:,2) = Cs_Xp_sp_sk;

Data_Xp_ps_map(:,:,3) = CNp_Xp_sp_sk;
Data_Xp_ps_map(:,:,4) = CNs_Xp_sp_sk;

Data_Xp_ps_map(:,:,5) = CNPp_Xp_sp_sk;
Data_Xp_ps_map(:,:,6) = CNPs_Xp_sp_sk;

Data_Xp_ps_map(301:360,:,:) = [];

% Net change in X
Cp_pol_sp_sk = nanmean(Cp_pol_maps(:,:,104:113),3) - nanmean(Cp_pol_maps(:,:,1:10),3);
Cs_pol_sp_sk = nanmean(Cs_pol_maps(:,:,104:113),3) - nanmean(Cs_pol_maps(:,:,1:10),3);

CNp_pol_sp_sk = nanmean(CNp_pol_maps(:,:,104:113),3) - nanmean(CNp_pol_maps(:,:,1:10),3);
CNs_pol_sp_sk = nanmean(CNs_pol_maps(:,:,104:113),3) - nanmean(CNs_pol_maps(:,:,1:10),3);

CNPp_pol_sp_sk = nanmean(CNPp_pol_maps(:,:,104:113),3) - nanmean(CNPp_pol_maps(:,:,1:10),3);
CNPs_pol_sp_sk = nanmean(CNPs_pol_maps(:,:,104:113),3) - nanmean(CNPs_pol_maps(:,:,1:10),3);

Data_X_ps_map(:,:,1) = Cp_pol_sp_sk;
Data_X_ps_map(:,:,2) = Cs_pol_sp_sk;

Data_X_ps_map(:,:,3) = CNp_pol_sp_sk;
Data_X_ps_map(:,:,4) = CNs_pol_sp_sk;

Data_X_ps_map(:,:,5) = CNPp_pol_sp_sk;
Data_X_ps_map(:,:,6) = CNPs_pol_sp_sk;

Data_X_ps_map(301:360,:,:) = [];

load('F:\My research\case2\JGR\Working\step1_reXcpx\map_X12.mat')
fig5 = figure
hold on
set(gcf,'position',[100 80 950 600])
axis off
maps_CNP = tight_subplot(3,2,[-0.08 -0.09],[0.06 0.02],[0 0.3])

Labels = {'(a)','(b)','(c)','(d)',...
          '(e)','(f)','(g)','(h)',...
          '(i)','(j)','(k)','(l)',...
          '(m)','(n)','(o)','(p)',...
          '(q)','(r)','(s)','(t)'};
      

% plant C sink
for i= 1:2:5
    i
    
    map_i = flipud(Data_Xp_ps_map(:,:,i))./1000;   % Unit: KgC m-2
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
    colormap(figCNP, map_X12)
    caxis([-5 5])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

% soil C sink
for i= 2:2:6
    i
    
    map_i = flipud(Data_Xp_ps_map(:,:,i))./1000;   % Unit: KgC m-2
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
    colormap(figCNP, map_X12)
    caxis([-5 5])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end
h = colorbar
h.Location = 'southoutside';
h.Position = [0.0572,0.087,0.59,0.0246]
h.Label.FontName = 'Arial'
h.Label.FontSize = 10;
text(-3.3712,-2.3437, '\DeltaXp (KgC m^-^2)','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',10)
    
text(-9.3102,10.4495, 'Plant:','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12)   
text(-2.5,10.4495, 'Soil:','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12)      
 
%% BarPlot for Net change in Xp and land C storage over 1901-2013
clearvars -except fig5
load('F:\My research\case2\JGR\1_Figures\Figure5\Plant_Soil_X_Xc_Xp.mat') 

% the spatial distribution of the magnitude of disequilibrium for plant and soil
% C-only
Cp_Xp_sp_sk = nanmean(Cp_Xp_maps(:,:,104:113),3) - nanmean(Cp_Xp_maps(:,:,1:10),3);
Cs_Xp_sp_sk = nanmean(Cs_Xp_maps(:,:,104:113),3) - nanmean(Cs_Xp_maps(:,:,1:10),3);
% CN
CNp_Xp_sp_sk = nanmean(CNp_Xp_maps(:,:,104:113),3) - nanmean(CNp_Xp_maps(:,:,1:10),3);
CNs_Xp_sp_sk = nanmean(CNs_Xp_maps(:,:,104:113),3) - nanmean(CNs_Xp_maps(:,:,1:10),3);
% CNP
CNPp_Xp_sp_sk = nanmean(CNPp_Xp_maps(:,:,104:113),3) - nanmean(CNPp_Xp_maps(:,:,1:10),3);
CNPs_Xp_sp_sk = nanmean(CNPs_Xp_maps(:,:,104:113),3) - nanmean(CNPs_Xp_maps(:,:,1:10),3);

Data_Xp_ps_map(:,:,1) = Cp_Xp_sp_sk;
Data_Xp_ps_map(:,:,2) = Cs_Xp_sp_sk;
Data_Xp_ps_map(:,:,3) = CNp_Xp_sp_sk;
Data_Xp_ps_map(:,:,4) = CNs_Xp_sp_sk;
Data_Xp_ps_map(:,:,5) = CNPp_Xp_sp_sk;
Data_Xp_ps_map(:,:,6) = CNPs_Xp_sp_sk;


% Net change in X
Cp_pol_sp_sk = nanmean(Cp_pol_maps(:,:,104:113),3) - nanmean(Cp_pol_maps(:,:,1:10),3);
Cs_pol_sp_sk = nanmean(Cs_pol_maps(:,:,104:113),3) - nanmean(Cs_pol_maps(:,:,1:10),3);
CNp_pol_sp_sk = nanmean(CNp_pol_maps(:,:,104:113),3) - nanmean(CNp_pol_maps(:,:,1:10),3);
CNs_pol_sp_sk = nanmean(CNs_pol_maps(:,:,104:113),3) - nanmean(CNs_pol_maps(:,:,1:10),3);
CNPp_pol_sp_sk = nanmean(CNPp_pol_maps(:,:,104:113),3) - nanmean(CNPp_pol_maps(:,:,1:10),3);
CNPs_pol_sp_sk = nanmean(CNPs_pol_maps(:,:,104:113),3) - nanmean(CNPs_pol_maps(:,:,1:10),3);

Data_X_ps_map(:,:,1) = Cp_pol_sp_sk;
Data_X_ps_map(:,:,2) = Cs_pol_sp_sk;
Data_X_ps_map(:,:,3) = CNp_pol_sp_sk;
Data_X_ps_map(:,:,4) = CNs_pol_sp_sk;
Data_X_ps_map(:,:,5) = CNPp_pol_sp_sk;
Data_X_ps_map(:,:,6) = CNPs_pol_sp_sk;

Data_reXp_PgG = [];
Data_reX_PgG = [];

for i = 1:6
    i
    reXp_gb = Data_Xp_ps_map(:,:,i).*cellarea .* 10^6./10^15;     % convert unit from gC m-2 into Pg C 
    reX_gb = Data_X_ps_map(:,:,i).*cellarea .* 10^6./10^15;
    
    Data_reXp_PgG(i) = nansum(reXp_gb(:));
    Data_reX_PgG(i) = nansum(reX_gb(:));
end



subplot_CNP = tight_subplot(2,1,[0.07 0.07],[0.12 0.07],[0.72 0.01])


Data_BarPlot_Xp_cnp = reshape(Data_reXp_PgG,[2,3])'
axes(subplot_CNP(1))
hold on
Xps_re = [1 3 5];

Bar_1P_2s = bar(Xps_re, Data_BarPlot_Xp_cnp, 'BarWidth',0.8)

Bar_1P_2s(1). FaceColor = [0.47,0.67,0.19];
Bar_1P_2s(2). FaceColor = [0.74,0.52,0.02];

Bar_1P_2s(1). EdgeColor = [0.24,0.41,0.02];
Bar_1P_2s(2). EdgeColor = [0.85,0.33,0.10];

Bar_1P_2s(1). LineWidth = 1.5;
Bar_1P_2s(2). LineWidth = 1.5;

plot([2 2],[0 450],'k--','LineWidth',1)
plot([4 4],[0 450],'k--','LineWidth',1)

set(gca,'linewidth',1.2,'box','on')
set(gca,'XLim',[0 6],'YLim',[0 450]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');

ylabel('\DeltaXp in plant/soil (PgC)','Fontname','Arial','FontSize',12)
xticks([1, 3, 5]);
xticklabels({'C-only','CN','CNP'})

legCNP = legend({'Plant','Soil'});    
set(legCNP,'color','none','EdgeColor','none','Fontname','Arial','Fontsize',9,'NumColumns',1,'Position',[0.91,0.861,0.083,0.0575])
text(0.4864,416.3, '(g)','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11)      

delete(subplot_CNP(2))

