clear
clc

% 1. Comparison of initial values, NPP0 and tuaE0, among C, CN and CNP
% 2. calculate changes of NPP and tuaE over 1901-2013
% 3. Comparison of three determinative components of Xc

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

C_NPP_maps = [];
C_tuaE_maps = [];

% read NPP data and convert into global map (unit: gC m-2 yr-1)
file_NPP = 'F:\My research\case2\working\cable-C\Matrix_out\';
% read ecosystem residence time data (unit: days)
file_tuaE = 'F:\My research\case2\working\cable-C\Transit_Matrix\';
ndays = 365;

for i = 1901:2013
    i
    
    C_NPPinf = xlsread([file_NPP, 'Carbon_1_',num2str(i),'_run_CN_cycle.csv']);  
    C_tuaEinf = xlsread([file_tuaE,'resident_All_',num2str(i),'.csv']);  
    
    C_NPP = C_NPPinf(:,1);
    C_tuaE = C_tuaEinf(:,1)./ndays;  % change unit from day to year
    
    % convert into global map based on mask
    mask_C_NPP = mask; index_NPP = 1;
    mask_C_tuaE = mask; index_tuaE = 1;
    
    for rowID=1:360;
        for colID=1:720
            
            if isnan(mask_C_NPP(rowID,colID)) == 0
               mask_C_NPP(rowID,colID) = C_NPP(index_NPP,1);
               index_NPP = index_NPP+1;
               index_NPP
            end
            
            if isnan(mask_C_tuaE(rowID,colID)) == 0
               mask_C_tuaE(rowID,colID) = C_tuaE(index_tuaE,1);
               index_tuaE = index_tuaE + 1;
            end

        end
    end
    
    C_NPP_maps(:,:,i-1900) = mask_C_NPP;
    C_tuaE_maps(:,:,i-1900) = mask_C_tuaE;
end

clearvars -except mask cellarea ...
                  C_NPP_maps C_tuaE_maps

%% CN

CN_NPP_maps = [];
CN_tuaE_maps = [];

% read NPP data and convert into global map (unit: gC m-2 yr-1)
file_NPP = 'F:\My research\case2\working\cable-CN\Matrix_out\';
% read ecosystem residence time data (unit: days)
file_tuaE = 'F:\My research\case2\working\cable-CN\Transit_Matrix\';
ndays = 365;

for i = 1901:2013
    i
    
    CN_NPPinf = xlsread([file_NPP, 'Carbon_1_',num2str(i),'_run_CN_cycle.csv']);  
    CN_tuaEinf = xlsread([file_tuaE,'resident_All_',num2str(i),'.csv']);  
    
    CN_NPP = CN_NPPinf(:,1);
    CN_tuaE = CN_tuaEinf(:,1)./ndays;  % change unit from day to year
    
    % convert into global map based on mask
    mask_CN_NPP = mask; index_NPP = 1;
    mask_CN_tuaE = mask; index_tuaE = 1;
    
    for rowID=1:360;
        for colID=1:720
            
            if isnan(mask_CN_NPP(rowID,colID)) == 0
               mask_CN_NPP(rowID,colID) = CN_NPP(index_NPP,1);
               index_NPP = index_NPP+1;
               index_NPP
            end
            
            if isnan(mask_CN_tuaE(rowID,colID)) == 0
               mask_CN_tuaE(rowID,colID) = CN_tuaE(index_tuaE,1);
               index_tuaE = index_tuaE + 1;
            end

        end
    end
    
    CN_NPP_maps(:,:,i-1900) = mask_CN_NPP;
    CN_tuaE_maps(:,:,i-1900) = mask_CN_tuaE;
end

clearvars -except mask cellarea ...
                  C_NPP_maps C_tuaE_maps ...
                  CN_NPP_maps CN_tuaE_maps

%% CNP

CNP_NPP_maps = [];
CNP_tuaE_maps = [];

% read NPP data and convert into global map (unit: gC m-2 yr-1)
file_NPP = 'F:\My research\case2\working\cable-CNP\Matrix_out\';
% read ecosystem residence time data (unit: days)
file_tuaE = 'F:\My research\case2\working\cable-CNP\Transit_Matrix\';
ndays = 365;

for i = 1901:2013
    i
    
    CNP_NPPinf = xlsread([file_NPP, 'Carbon_1_',num2str(i),'_run_CN_cycle.csv']);  
    CNP_tuaEinf = xlsread([file_tuaE,'resident_All_',num2str(i),'.csv']);  
    
    CNP_NPP = CNP_NPPinf(:,1);
    CNP_tuaE = CNP_tuaEinf(:,1)./ndays;  % change unit from day to year
    
    % convert into global map based on mask
    mask_CNP_NPP = mask; index_NPP = 1;
    mask_CNP_tuaE = mask; index_tuaE = 1;
    
    for rowID=1:360;
        for colID=1:720
            
            if isnan(mask_CNP_NPP(rowID,colID)) == 0
               mask_CNP_NPP(rowID,colID) = CNP_NPP(index_NPP,1);
               index_NPP = index_NPP+1;
               index_NPP
            end
            
            if isnan(mask_CNP_tuaE(rowID,colID)) == 0
               mask_CNP_tuaE(rowID,colID) = CNP_tuaE(index_tuaE,1);
               index_tuaE = index_tuaE + 1;
            end

        end
    end
    
    CNP_NPP_maps(:,:,i-1900) = mask_CNP_NPP;
    CNP_tuaE_maps(:,:,i-1900) = mask_CNP_tuaE;
end

clearvars -except mask cellarea ...
                  C_NPP_maps C_tuaE_maps ...
                  CN_NPP_maps CN_tuaE_maps ...
                  CNP_NPP_maps CNP_tuaE_maps 
              
cd('F:\My research\case2\JGR\Working\step3_reNPPtuaE_fct')      
save 'step3_reNPPtuaE_fct_data.mat'
load ('F:\My research\case2\JGR\Working\step3_reNPPtuaE_fct\step3_reNPPtuaE_fct_data.mat')              
%% FigureS
% initial values for NPP and Ecosystem residence time
% estimated as the initial mean of first decade.

% C-only
% NPP
C_NPP_10bg = C_NPP_maps(:,:,1:10);
C_NPP_10ed = C_NPP_maps(:,:,104:113);
C_NPP_10bg_ag = nanmean(C_NPP_10bg,3);
C_NPP_10ed_ag = nanmean(C_NPP_10ed,3);

Net_C_NPP_map = C_NPP_10ed_ag - C_NPP_10bg_ag;

% tuaE
C_tuaE_10bg = C_tuaE_maps(:,:,1:10);
C_tuaE_10ed = C_tuaE_maps(:,:,104:113);
C_tuaE_10bg_ag = nanmean(C_tuaE_10bg,3);
C_tuaE_10ed_ag = nanmean(C_tuaE_10ed,3);

Net_C_tuaE_map = C_tuaE_10ed_ag - C_tuaE_10bg_ag;

% CN
% NPP
CN_NPP_10bg = CN_NPP_maps(:,:,1:10);
CN_NPP_10ed = CN_NPP_maps(:,:,104:113);
CN_NPP_10bg_ag = nanmean(CN_NPP_10bg,3);
CN_NPP_10ed_ag = nanmean(CN_NPP_10ed,3);

Net_CN_NPP_map = CN_NPP_10ed_ag - CN_NPP_10bg_ag;

% tuaE
CN_tuaE_10bg = CN_tuaE_maps(:,:,1:10);
CN_tuaE_10ed = CN_tuaE_maps(:,:,104:113);
CN_tuaE_10bg_ag = nanmean(CN_tuaE_10bg,3);
CN_tuaE_10ed_ag = nanmean(CN_tuaE_10ed,3);

Net_CN_tuaE_map = CN_tuaE_10ed_ag - CN_tuaE_10bg_ag;

% CNP
% NPP
CNP_NPP_10bg = CNP_NPP_maps(:,:,1:10);
CNP_NPP_10ed = CNP_NPP_maps(:,:,104:113);
CNP_NPP_10bg_ag = nanmean(CNP_NPP_10bg,3);
CNP_NPP_10ed_ag = nanmean(CNP_NPP_10ed,3);

Net_CNP_NPP_map = CNP_NPP_10ed_ag - CNP_NPP_10bg_ag;
 
% tuaE
CNP_tuaE_10bg = CNP_tuaE_maps(:,:,1:10);
CNP_tuaE_10ed = CNP_tuaE_maps(:,:,104:113);
CNP_tuaE_10bg_ag = nanmean(CNP_tuaE_10bg,3);
CNP_tuaE_10ed_ag = nanmean(CNP_tuaE_10ed,3);

Net_CNP_tuaE_map = CNP_tuaE_10ed_ag - CNP_tuaE_10bg_ag;


% preparing for data 
DataMap_NPPtuaE(:,:,1) = C_NPP_10bg_ag;
DataMap_NPPtuaE(:,:,2) = C_tuaE_10bg_ag;

DataMap_NPPtuaE(:,:,3) = CN_NPP_10bg_ag;
DataMap_NPPtuaE(:,:,4) = CN_tuaE_10bg_ag;

DataMap_NPPtuaE(:,:,5) = CNP_NPP_10bg_ag;
DataMap_NPPtuaE(:,:,6) = CNP_tuaE_10bg_ag;

DataPoint_NPPtuaE0 = DataMap_NPPtuaE;

DataMap_NPPtuaE(301:360,:,:) = [];


%ax = gca;
%mymap_CNP_tuaE2 = colormap(ax)
%mymap_CNP_tuaE2 = mymap_CNP_tuaE2(1:16,:);
%save('mymap_CNP_tuaE2','mymap_CNP_tuaE2')

% global maps for initial state of NPP and tuaE
load('F:\My research\case2\JGR\Working\step4_PlantSoil\inCtuaE0_soilplant\mymap_inC0_sp.mat')
load('F:\My research\case2\JGR\Working\step2_NPP0_tuaE0\mymap_CNP_tuaE2.mat')
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
    
    map_i = flipud(DataMap_NPPtuaE(:,:,i))./1000;   % Unit: KgC m-2 yr-1
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

h = colorbar
h.Label.String = 'Initial NPP (KgC m^-^2 yr^-^1)'
h.Label.FontName = 'Arial'
h.Label.FontSize = 12;


for i= 2:2:6
    i
    
    map_i = flipud(DataMap_NPPtuaE(:,:,i));   % Unit: year
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
    colormap(figCNP,mymap_CNP_tuaE2)
    caxis([0 800])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

h = colorbar
h.Label.String = 'Initial ecosystem residence time (year)'
h.Label.FontName = 'Arial'
h.Label.FontSize = 12;

% global plot for cells 

for i = 1:2:5
    i
    
    DataPoint_NPPtuaE0(:,:,i) = DataPoint_NPPtuaE0(:,:,i)./1000;  % convert unit of NPP into KgC m-2
    
end

% NPP2tuaE 
C_NPP0 = DataPoint_NPPtuaE0(:,:,1);
C_tuaE0 = DataPoint_NPPtuaE0(:,:,2);

CN_NPP0 = DataPoint_NPPtuaE0(:,:,3);
CN_tuaE0 = DataPoint_NPPtuaE0(:,:,4);

CNP_NPP0 = DataPoint_NPPtuaE0(:,:,5);
CNP_tuaE0 = DataPoint_NPPtuaE0(:,:,6);

figure
hold on
plot(C_NPP0(:),C_tuaE0(:),...
     'Marker','.','MarkerEdgeColor',[0.34 0.33 0.33],...
     'MarkerFaceColor',[0.34 0.33 0.33],...
     'MarkerSize', 6,...
     'LineStyle','none')
alpha(0.9)

plot(CN_NPP0(:),CN_tuaE0(:),...
     'Marker','.','MarkerEdgeColor',[1 0 0],...
     'MarkerFaceColor',[1 0 0],...
     'MarkerSize', 6,...
     'LineStyle','none')
alpha(0.5)
 
plot(CNP_NPP0(:),CNP_tuaE0(:),...
     'Marker','.','MarkerEdgeColor',[0 0 1],...
     'MarkerFaceColor',[0 0 1],...
     'MarkerSize', 4,...
     'LineStyle','none')
alpha(0.25)

set(gca,'linewidth',1.2,'box','on')
set(gca,'XLim',[0 3],'YLim',[0 1200])
set(gca,'Fontname','Arial','FontSize',12);

ylabel('Ecosystem residence time (year)','Fontname','Arial','FontSize',14)
xlabel('NPP (KgC m^-^1)','Fontname','Arial','FontSize',14)

text(0.15,1100,'(a)','Fontname','Arial','FontSize',14)

figure
hold on

plot(nanmean(C_NPP0(:)),nanmean(C_tuaE0(:)),...
     'Marker','v','MarkerEdgeColor',[0.34 0.33 0.33],...
     'MarkerFaceColor',[0.34 0.33 0.33],...
     'MarkerSize', 14,...
     'LineStyle','none')


plot(nanmean(CN_NPP0(:)),nanmean(CN_tuaE0(:)),...
     'Marker','o','MarkerEdgeColor',[1 0 0],...
     'MarkerFaceColor',[1 0 0],...
     'MarkerSize', 14,...
     'LineStyle','none')
kCN =  nanmean(CN_NPP0(:))* nanmean(CN_tuaE0(:));

plot(nanmean(CNP_NPP0(:)),nanmean(CNP_tuaE0(:)),...
     'Marker','s','MarkerEdgeColor',[0 0 1],...
     'MarkerFaceColor',[0 0 1],...
     'MarkerSize', 14,...
     'LineStyle','none') 
kCNP =  nanmean(CNP_NPP0(:))* nanmean(CNP_tuaE0(:)); 
 
 
set(gca,'linewidth',1.2,'box','off')
set(gca,'XLim',[0.2 0.5],'YLim',[95 102])
set(gca,'Fontname','Arial','FontSize',12); 

text(nanmean(C_NPP0(:))-0.03,nanmean(C_tuaE0(:))+1,'C-only','Fontname','Arial','FontSize',11)
text(nanmean(CN_NPP0(:))+0.015,nanmean(CN_tuaE0(:))+0.5,'CN','Fontname','Arial','FontSize',11)
text(nanmean(CNP_NPP0(:))-0.03,nanmean(CNP_tuaE0(:))-1,'CNP','Fontname','Arial','FontSize',11)

text(0.23,101.3,'(b)','Fontname','Arial','FontSize',14) 
%% FigureS
% Nutrient limitation on initial NPP and tuaE

Nlimt_NPP0 = CN_NPP_10bg_ag - C_NPP_10bg_ag;
Nlimt_tuaE0 = CN_tuaE_10bg_ag - C_tuaE_10bg_ag;

Plimt_NPP0 = CNP_NPP_10bg_ag - CN_NPP_10bg_ag;
Plimt_tuaE0 = CNP_tuaE_10bg_ag - CN_tuaE_10bg_ag;

% preparing for data
Data_NPlimt0(:,:,1) = Nlimt_NPP0;
Data_NPlimt0(:,:,2) = Nlimt_tuaE0;

Data_NPlimt0(:,:,3) = Plimt_NPP0;
Data_NPlimt0(:,:,4) = Plimt_tuaE0;

Data_NPlimt0(301:360,:,:) = [];

%ax = gca;
%mymap_NP_limt0 = colormap(ax)
%save('mymap_NP_limt0','mymap_NP_limt0')

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
    
    map_i = flipud(Data_NPlimt0(:,:,i))./1000;   % Unit: KgC m-2 yr-1
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
    
    map_i = flipud(Data_NPlimt0(:,:,i));   % Unit: years
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
    caxis([-30 30])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end      
  
     
%% Figure3
% re_NPP * tuaE0, re_tuaE * NPP0,  re_NPP * re_tuaE

% C-only
fct_reNPP_C = Net_C_NPP_map.* C_tuaE_10bg_ag;
fct_retuaE_C = Net_C_tuaE_map.* C_NPP_10bg_ag;
fct_reNPPtuaE_C = Net_C_NPP_map.* Net_C_tuaE_map;

% CN
fct_reNPP_CN = Net_CN_NPP_map.* CN_tuaE_10bg_ag;
fct_retuaE_CN = Net_CN_tuaE_map.* CN_NPP_10bg_ag;
fct_reNPPtuaE_CN = Net_CN_NPP_map.* Net_CN_tuaE_map;

% CNP
fct_reNPP_CNP = Net_CNP_NPP_map.* CNP_tuaE_10bg_ag;
fct_retuaE_CNP = Net_CNP_tuaE_map.* CNP_NPP_10bg_ag;
fct_reNPPtuaE_CNP = Net_CNP_NPP_map.* Net_CNP_tuaE_map;

% preparing for dataset 
Data_reNPPtuaE(:,:,1) = fct_reNPP_C;
Data_reNPPtuaE(:,:,2) = fct_retuaE_C;
Data_reNPPtuaE(:,:,3) = fct_reNPPtuaE_C;

Data_reNPPtuaE(:,:,4) = fct_reNPP_CN;
Data_reNPPtuaE(:,:,5) = fct_retuaE_CN;
Data_reNPPtuaE(:,:,6) = fct_reNPPtuaE_CN;

Data_reNPPtuaE(:,:,7) = fct_reNPP_CNP;
Data_reNPPtuaE(:,:,8) = fct_retuaE_CNP;
Data_reNPPtuaE(:,:,9) = fct_reNPPtuaE_CNP;

Data_reNPPtuaE(301:360,:,:) = [];


load('F:\My research\case2\JGR\Working\step3_reNPPtuaE_fct\mymap_CNP_re.mat')
figure
set(gcf,'position',[100 100 750 450])
maps_CNP = tight_subplot(3,3,[-0.08 -0.08],[0.06 0.02],[0.01 0.01])
Labels = {'(a)','(b)','(c)','(d)',...
          '(e)','(f)','(g)','(h)',...
          '(i)','(j)','(k)','(l)',...
          '(m)','(n)','(o)','(p)',...
          '(q)','(r)','(s)','(t)'};     
for i= 1:3:7
    i
    
    map_i = flipud(Data_reNPPtuaE(:,:,i))./1000;   % Unit: KgC m-2 
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
h = colorbar
h.Location = 'southoutside';
h.Position = [0.0624,0.0843,0.275,0.0247]
h.Label.FontName = 'Arial'
h.Label.FontSize = 12;
text(-0.0379,-2.39, 'KgC m^-^2','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',10)


for i= 2:3:8
    i
    
    map_i = flipud(Data_reNPPtuaE(:,:,i))./1000;   % Unit: KgC m-2
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


for i= 3:3:9
    i
    
    map_i = flipud(Data_reNPPtuaE(:,:,i))./1000;   % Unit: KgC m-2 yr-1
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
h = colorbar
h.Location = 'southoutside';
h.Position = [0.3837,0.0843,0.5243,0.0247]
h.Label.FontName = 'Arial'
h.Label.FontSize = 10;
text(-3.629,-2.39, 'KgC m^-^2','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',10)


text(-14.0432,10.4495, '\tau_E_0 X \DeltaNPP','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',10)        
text(-6.78,10.4495, '\Delta\tau_E X NPP_0','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',10)    
text(0.2,10.4495, '\Delta\tau_E X \DeltaNPP','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',10) 
    
text(-18,9, 'C','HorizontalAlignment','left',...
        'FontName','Arial','FontSize',10)        
text(-18,5, 'CN','HorizontalAlignment','left',...
        'FontName','Arial','FontSize',10)    
text(-18,1, 'CNP','HorizontalAlignment','left',...
        'FontName','Arial','FontSize',10)     
    
    
    


%% FigrueS of reNPP and retuaE under three different coupling schemes

% preparing for data 
DataMap_NPPtuaE_re(:,:,1) = Net_C_NPP_map;
DataMap_NPPtuaE_re(:,:,2) = Net_C_tuaE_map;

DataMap_NPPtuaE_re(:,:,3) = Net_CN_NPP_map;
DataMap_NPPtuaE_re(:,:,4) = Net_CN_tuaE_map;

DataMap_NPPtuaE_re(:,:,5) = Net_CNP_NPP_map;
DataMap_NPPtuaE_re(:,:,6) = Net_CNP_tuaE_map;

DataMap_NPPtuaE_re(301:360,:,:) = [];

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
    
    map_i = flipud(DataMap_NPPtuaE_re(:,:,i))./1000;   % Unit: KgC m-2 yr-1
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
    caxis([-0.5 0.5])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

for i= 2:2:6
    i
    
    map_i = flipud(DataMap_NPPtuaE_re(:,:,i));   % Unit: year
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
    colormap(figCNP,mymap_CNP_re)
    caxis([-3 3])
    set(gca,'box','off')
    axis off
    colorbar('off')

    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end







