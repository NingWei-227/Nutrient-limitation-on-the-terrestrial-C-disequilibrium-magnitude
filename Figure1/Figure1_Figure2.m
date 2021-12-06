clear
clc

% codes for Figure1 and Figure2
% Comparison of Xp, Xc and X under C-only, C-N and C-N-P interations
% read mask and cellarea data used to convert model outputs into global map
[num,text,raw]=xlsread('D:\Data_WN\OneDrive\W_cases\case2\case2\working\ALL\codes\yanmo0.5.csv');
for rowID=1:360
    for colID =1:720
        if raw{rowID,colID}=='NA'
            raw{rowID,colID}=-9999;
        end  
    end
end
mask = cell2mat(raw);
mask(mask==-9999) =  nan;
cellarea = xlsread('D:\Data_WN\OneDrive\W_cases\case2\case2\working\ALL\codes\area0.5.csv');   % unit: km2 = 10^6 m2
clearvars -except mask cellarea

%% C-only 
cd('D:\Data_WN\OneDrive\W_cases\case2\case2\working\cable-C\Transit_Matrix');
C_cap_maps = [];    % define a null matrix to store data of carbon storage capacity
C_pot_maps = [];    % define a null matrix to store data of disequilibrium magnitude
C_pol_maps = [];    % define a null matrix to store data of carbon pool

% read data of C storege capacity, C storage and disequilibrium and convert into global maps (unit: )
for i = 1901:2013
    i
    
    C_cap9 = xlsread(['C_capacity_',num2str(i),'.csv']);   % C storage capacity of 9 C pools (unit: gC m-2)
    C_pot9 = xlsread(['C_potencial_',num2str(i),'.csv']);  % the disequilibrium 
    C_pol9 = xlsread(['X_xx_',num2str(i),'.csv']);         % C storage of 9 C pools
    
    % sum up 9 C pools
    C_capE = sum(C_cap9,2);                               
    C_potE = sum(C_pot9,2);
    C_polE = sum(C_pol9,2);
    
    % convert into global map based on mask
    mask_C_cap = mask; index_cap = 1;
    mask_C_pol = mask; index_pol = 1;
    mask_C_pot = mask; index_pot = 1;
    
    for rowID=1:360;
        for colID=1:720
            
            if isnan(mask_C_cap(rowID,colID)) == 0
               mask_C_cap(rowID,colID) = C_capE(index_cap,1);
               index_cap = index_cap+1;
               index_cap
            end
            if isnan(mask_C_pol(rowID,colID)) == 0
               mask_C_pol(rowID,colID) = C_polE(index_pol,1);
               index_pol = index_pol + 1;
            end
            if isnan(mask_C_pot(rowID,colID)) == 0
               mask_C_pot(rowID,colID) = C_potE(index_pot,1);
               index_pot = index_pot +1;
            end
            
        end
    end
    
    C_cap_maps(:,:,i-1900) = mask_C_cap;
    C_pot_maps(:,:,i-1900) = mask_C_pot;
    C_pol_maps(:,:,i-1900) = mask_C_pol;
      
end

% global estimation of Cstorage capacity, Cstorage and disequilibrium (unit: Pg C)
C_cap_yr = [];
C_pot_yr = [];
C_pol_yr = [];

for yr = 1:113
    yr 
    
    C_cap_gb=  C_cap_maps(:,:,yr).*cellarea .* 10^6./10^15;    % convert unit from gC m-2 into Pg C
    C_pot_gb = C_pot_maps(:,:,yr).*cellarea .* 10^6./10^15;
    C_pol_gb = C_pol_maps(:,:,yr).*cellarea .* 10^6./10^15;
    
    C_cap_yr(yr)  = nansum(C_cap_gb(:));
    C_pot_yr(yr)  = nansum(C_pot_gb(:));
    C_pol_yr(yr)  = nansum(C_pol_gb(:));
    
end


clearvars -except C_cap_maps C_pol_maps C_pot_maps ...
                  C_cap_yr C_pot_yr C_pol_yr ...
                  mask cellarea

%% CABLE-CN
cd('D:\Data_WN\OneDrive\W_cases\case2\case2\working\cable-CN\Transit_Matrix');

CN_cap_maps = [];
CN_pot_maps = [];
CN_pol_maps = [];

% read data of C storege capacity, C storage and disequilibrium and convert
% into global maps (unit: )
for i = 1901:2013
    i
    
    CN_cap9 = xlsread(['C_capacity_',num2str(i),'.csv']);   % C storage capacity of 9 C pools (unit: gC m-2)
    CN_pot9 = xlsread(['C_potencial_',num2str(i),'.csv']);  % the disequilibrium 
    CN_pol9 = xlsread(['X_xx_',num2str(i),'.csv']);         % C storage of 9 C pools
    
    % sum up 9 C pools
    CN_capE = sum(CN_cap9,2);                               
    CN_potE = sum(CN_pot9,2);
    CN_polE = sum(CN_pol9,2);
    
    % convert into global map based on mask
    mask_CN_cap = mask; index_cap = 1;
    mask_CN_pol = mask; index_pol = 1;
    mask_CN_pot = mask; index_pot = 1;
    
    for rowID=1:360
        for colID=1:720
            
            if isnan(mask_CN_cap(rowID,colID)) == 0
               mask_CN_cap(rowID,colID) = CN_capE(index_cap,1);
               index_cap = index_cap+1;
               index_cap
            end
            if isnan(mask_CN_pol(rowID,colID)) == 0 
               mask_CN_pol(rowID,colID) = CN_polE(index_pol,1);
               index_pol = index_pol + 1;
            end
            if isnan(mask_CN_pot(rowID,colID)) == 0
               mask_CN_pot(rowID,colID) = CN_potE(index_pot,1);
               index_pot = index_pot +1;
            end
            
        end
    end
    
    CN_cap_maps(:,:,i-1900) = mask_CN_cap;
    CN_pot_maps(:,:,i-1900) = mask_CN_pot;
    CN_pol_maps(:,:,i-1900) = mask_CN_pol;
      
end


% global estimation of Cstorage capacity, Cstorage and disequilibrium (unit: Pg C)
CN_cap_yr = [];
CN_pot_yr = [];
CN_pol_yr = [];
for yr = 1:113
    yr 
    
    CN_cap_gb=  CN_cap_maps(:,:,yr).*cellarea .* 10^6./10^15;     % convert unit from gC m-2 into Pg C
    CN_pot_gb = CN_pot_maps(:,:,yr).*cellarea .* 10^6./10^15;
    CN_pol_gb = CN_pol_maps(:,:,yr).*cellarea .* 10^6./10^15;
    
    CN_cap_yr(yr)  = nansum(CN_cap_gb(:));
    CN_pot_yr(yr)  = nansum(CN_pot_gb(:));
    CN_pol_yr(yr)  = nansum(CN_pol_gb(:));
    
end

clearvars -except C_cap_maps C_pol_maps C_pot_maps ...
                  C_cap_yr C_pot_yr C_pol_yr ...
                  mask cellarea ...
                  CN_cap_maps CN_pol_maps CN_pot_maps ...
                  CN_cap_yr CN_pot_yr CN_pol_yr
                 
%% CABLE-CNP
cd('D:\Data_WN\OneDrive\W_cases\case2\case2\working\cable-CNP\Transit_Matrix');

CNP_cap_maps = [];
CNP_pot_maps = [];
CNP_pol_maps = [];

% read data of C storege capacity, C storage and disequilibrium and convert
% into global maps (unit: )
for i = 1901:2013
    i
    
    CNP_cap9 = xlsread(['C_capacity_',num2str(i),'.csv']);   % C storage capacity of 9 C pools (unit: gC m-2)
    CNP_pot9 = xlsread(['C_potencial_',num2str(i),'.csv']);  % the disequilibrium 
    CNP_pol9 = xlsread(['X_xx_',num2str(i),'.csv']);         % C storage of 9 C pools
    
    % sum up 9 C pools
    CNP_capE = sum(CNP_cap9,2);                               
    CNP_potE = sum(CNP_pot9,2);
    CNP_polE = sum(CNP_pol9,2);
    
    % convert into global map based on mask
    mask_CNP_cap = mask; index_cap = 1;
    mask_CNP_pol = mask; index_pol = 1;
    mask_CNP_pot = mask; index_pot = 1;
    
    for rowID=1:360
        for colID=1:720
            
            if isnan(mask_CNP_cap(rowID,colID)) == 0
               mask_CNP_cap(rowID,colID) = CNP_capE(index_cap,1);
               index_cap = index_cap+1;
               index_cap
            end
            if isnan(mask_CNP_pol(rowID,colID)) == 0
               mask_CNP_pol(rowID,colID) = CNP_polE(index_pol,1);
               index_pol = index_pol + 1;
            end
            if isnan(mask_CNP_pot(rowID,colID)) == 0
               mask_CNP_pot(rowID,colID) = CNP_potE(index_pot,1);
               index_pot = index_pot +1;
            end
            
        end
    end
    
    CNP_cap_maps(:,:,i-1900) = mask_CNP_cap;
    CNP_pot_maps(:,:,i-1900) = mask_CNP_pot;
    CNP_pol_maps(:,:,i-1900) = mask_CNP_pol;
      
end


% global estimation of Cstorage capacity, Cstorage and disequilibrium (unit: Pg C)
CNP_cap_yr = [];
CNP_pot_yr = [];
CNP_pol_yr = [];
for yr = 1:113
    yr 
    
    CNP_cap_gb=  CNP_cap_maps(:,:,yr).*cellarea .* 10^6./10^15;     % convert unit from gC m-2 into Pg C
    CNP_pot_gb = CNP_pot_maps(:,:,yr).*cellarea .* 10^6./10^15;     % cellarea umnit: km2
    CNP_pol_gb = CNP_pol_maps(:,:,yr).*cellarea .* 10^6./10^15;
    
    CNP_cap_yr(yr)  = nansum(CNP_cap_gb(:));
    CNP_pot_yr(yr)  = nansum(CNP_pot_gb(:));
    CNP_pol_yr(yr)  = nansum(CNP_pol_gb(:));
    
end

clearvars -except C_cap_maps C_pol_maps C_pot_maps ...
                  C_cap_yr C_pot_yr C_pol_yr ...
                  mask cellarea ...
                  CN_cap_maps CN_pol_maps CN_pot_maps ...
                  CN_cap_yr CN_pot_yr CN_pol_yr ...
                  CNP_cap_maps CNP_pol_maps CNP_pot_maps ...
                  CNP_cap_yr CNP_pot_yr CNP_pol_yr
              
save D:\Data_WN\OneDrive\W_cases\case2\case2\writing\JGR\Resubmit\1_Figure\Figure1\Figure1mataData.mat
load D:\Data_WN\OneDrive\W_cases\case2\case2\writing\JGR\Resubmit\1_Figure\Figure1\Figure1mataData.mat
%% Figure 1 
figure      
set(gcf,'position',[100 100 970 420])
CNP_fig = tight_subplot(1,2,[0 0.08],[0.15 0.02],[0.1 0.02])

% panel (a): Temporal dynamics of Xc, X and Xp
axes(CNP_fig(1));
hold on
X_cnp(:,1) = C_pol_yr;
X_cnp(:,2) = CN_pol_yr;
X_cnp(:,3) = CNP_pol_yr;

lower(1:113,:) = 0;
shade(:,:,1) = [lower C_pot_yr'];
shade(:,:,2) = [lower CN_pot_yr'];
shade(:,:,3) = [lower CNP_pot_yr'];

colorCNP = [0 0 0; 1 0 0; 0 0 1];
Years = 1901:2013;

CNP_leg = boundedline(Years',X_cnp,shade,'cmap',colorCNP,'alpha');
set(CNP_leg,'LineWidth',1.8);
set(gca,'linewidth',1.2,'box','on')
set(gca,'XLim',[1900 2013],'YLim',[1300 4200])
set(gca,'Fontname','Arial','FontSize',12);
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');

legCNP = legend(CNP_leg,{'C-only', 'CN', 'CNP'})
set(legCNP,'color','none','EdgeColor','none','Fontname','Arial','Fontsize',12,'NumColumns',1)
ylabel('Land C storage (PgC)','Fontname','Arial','FontSize',14)
%xlabel('Years','Fontname','Arial','FontSize',14)
text(1905, 4000,'(a)','Fontname','Arial','FontSize',14)

% panel (b): absolute change in Cstorage, Cstorage capacity and disequilibrium    
% Units: Pg C [bars]
C_capYr_ChgRe = nanmean(C_cap_yr(104:113)) - nanmean(C_cap_yr(1:10));
C_potYr_ChgRe = nanmean(C_pot_yr(104:113)) - nanmean(C_pot_yr(1:10));
C_polYr_ChgRe = nanmean(C_pol_yr(104:113)) - nanmean(C_pol_yr(1:10));

CN_capYr_ChgRe = nanmean(CN_cap_yr(104:113)) - nanmean(CN_cap_yr(1:10));
CN_potYr_ChgRe = nanmean(CN_pot_yr(104:113)) - nanmean(CN_pot_yr(1:10));
CN_polYr_ChgRe = nanmean(CN_pol_yr(104:113)) - nanmean(CN_pol_yr(1:10));

CNP_capYr_ChgRe = nanmean(CNP_cap_yr(104:113)) - nanmean(CNP_cap_yr(1:10));
CNP_potYr_ChgRe = nanmean(CNP_pot_yr(104:113)) - nanmean(CNP_pot_yr(1:10));
CNP_polYr_ChgRe = nanmean(CNP_pol_yr(104:113)) - nanmean(CNP_pol_yr(1:10));

axes(CNP_fig(2));
hold on
leg(1) = bar([1 5 9], [C_potYr_ChgRe C_capYr_ChgRe C_polYr_ChgRe],'barWidth',0.2,...
    'FaceColor',[0.34 0.33 0.33],'EdgeColor',[0.34 0.33 0.33]);
leg(2) = bar([2 6 10], [CN_potYr_ChgRe CN_capYr_ChgRe CN_polYr_ChgRe],'barWidth',0.2,...
    'FaceColor',[1 0 0],'EdgeColor',[1 0 0]);
leg(3) = bar([3 7 11], [CNP_potYr_ChgRe CNP_capYr_ChgRe CNP_polYr_ChgRe],'barWidth',0.2,...
    'FaceColor',[0 0 1],'EdgeColor',[0 0 1]);
plot([4 4],[0 800],'k--','LineWidth',1)
plot([8 8],[0 800],'k--','LineWidth',1)

set(gca,'linewidth',1.2,'box','on')
set(gca,'XLim',[0 12],'YLim',[0 800])
set(gca,'Fontname','Arial','FontSize',12);
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');

legCNP = legend(leg,{'C','CN','CNP'});    
set(legCNP,'color','none','EdgeColor','none','Fontname','Arial','Fontsize',10,'NumColumns',1)
ylabel('Absolute changes of X_p, X_c and X (PgC)','Fontname','Arial','FontSize',14)
xticks([2 6 10]);
xticklabels({'\DeltaXp','\DeltaXc','\DeltaX'})
text(0.5, 750,'(b)','Fontname','Arial','FontSize',14)


%% Figure2
% Net change in the magnitude of disequilibrium under three different coupling schemes 
% Units for spatial data: KgC m-2
% the difference between the last 10 years and initial 10 years

% C-only
C_pol_10bg = C_pol_maps(:,:,1:10);       % 1901-1910 Cstorage data
C_pol_10ed = C_pol_maps(:,:,104:113);    % 2004-2013 Cstorage data
C_cap_10bg = C_cap_maps(:,:,1:10);       % 1901-1910 Cstorage capacity
C_cap_10ed = C_cap_maps(:,:,104:113);    % 2004-2013 Cstorage capacity
C_pot_10bg = C_pot_maps(:,:,1:10);       % 1901-1910 disequilibrium
C_pot_10ed = C_pot_maps(:,:,104:113);    % 2004-2013 disequilibrium
% decadal mean
C_pol_10bg_ag = nanmean(C_pol_10bg,3);
C_pol_10ed_ag = nanmean(C_pol_10ed,3);
C_cap_10bg_ag = nanmean(C_cap_10bg,3);
C_cap_10ed_ag = nanmean(C_cap_10ed,3);
C_pot_10bg_ag = nanmean(C_pot_10bg,3);
C_pot_10ed_ag = nanmean(C_pot_10ed,3);
% Absolute change in Cstorage, Cstorage capacity and disequilibrium
Net_C_pol = C_pol_10ed_ag - C_pol_10bg_ag;
Net_C_cap = C_cap_10ed_ag - C_cap_10bg_ag;
Net_C_pot = C_pot_10ed_ag - C_pot_10bg_ag;

% CN
CN_pol_10bg = CN_pol_maps(:,:,1:10);      % 1901-1910 Cstorage data
CN_pol_10ed = CN_pol_maps(:,:,104:113);   % 2004-2013 Cstorage data 
CN_cap_10bg = CN_cap_maps(:,:,1:10);      % 1901-1910 Cstorage capacity     
CN_cap_10ed = CN_cap_maps(:,:,104:113);   % 2004-2013 Cstorage capacity
CN_pot_10bg = CN_pot_maps(:,:,1:10);      % 1901-1910 disequilibrium
CN_pot_10ed = CN_pot_maps(:,:,104:113);   % 2004-2013 disequilibrium
% decadal mean
CN_pol_10bg_ag = nanmean(CN_pol_10bg,3);
CN_pol_10ed_ag = nanmean(CN_pol_10ed,3);
CN_cap_10bg_ag = nanmean(CN_cap_10bg,3);
CN_cap_10ed_ag = nanmean(CN_cap_10ed,3);
CN_pot_10bg_ag = nanmean(CN_pot_10bg,3);
CN_pot_10ed_ag = nanmean(CN_pot_10ed,3);
% Absolute change in Cstorage, Cstorage capacity and disequilibrium
Net_CN_pol = CN_pol_10ed_ag - CN_pol_10bg_ag;
Net_CN_cap = CN_cap_10ed_ag - CN_cap_10bg_ag;
Net_CN_pot = CN_pot_10ed_ag - CN_pot_10bg_ag;

% CNP
CNP_pol_10bg = CNP_pol_maps(:,:,1:10);        % 1901-1910 Cstorage data
CNP_pol_10ed = CNP_pol_maps(:,:,104:113);     % 2004-2013 Cstorage data 
CNP_cap_10bg = CNP_cap_maps(:,:,1:10);        % 1901-1910 Cstorage capacity  
CNP_cap_10ed = CNP_cap_maps(:,:,104:113);     % 2004-2013 Cstorage capacity
CNP_pot_10bg = CNP_pot_maps(:,:,1:10);        % 1901-1910 disequilibrium
CNP_pot_10ed = CNP_pot_maps(:,:,104:113);     % 2004-2013 disequilibrium
% decadal mean
CNP_pol_10bg_ag = nanmean(CNP_pol_10bg,3);
CNP_pol_10ed_ag = nanmean(CNP_pol_10ed,3);
CNP_cap_10bg_ag = nanmean(CNP_cap_10bg,3);
CNP_cap_10ed_ag = nanmean(CNP_cap_10ed,3);
CNP_pot_10bg_ag = nanmean(CNP_pot_10bg,3);
CNP_pot_10ed_ag = nanmean(CNP_pot_10ed,3);
% Absolute change in Cstorage, Cstorage capacity and disequilibrium
Net_CNP_pol = CNP_pol_10ed_ag - CNP_pol_10bg_ag;
Net_CNP_cap = CNP_cap_10ed_ag - CNP_cap_10bg_ag;
Net_CNP_pot = CNP_pot_10ed_ag - CNP_pot_10bg_ag;

Data_CNP_map(:,:,1) = Net_C_pot;
Data_CNP_map(:,:,2) = Net_CN_pot;
Data_CNP_map(:,:,3) = Net_CNP_pot;

Data_CNP_map(:,:,4) = Net_C_cap;
Data_CNP_map(:,:,5) = Net_CN_cap;
Data_CNP_map(:,:,6) = Net_CNP_cap;

Data_CNP_map(:,:,7) = Net_C_pol;
Data_CNP_map(:,:,8) = Net_CN_pol;
Data_CNP_map(:,:,9) = Net_CNP_pol;

Data_CNP_map = Data_CNP_map./1000;  % Unit: KgC m-2
Data_CNP_map(301:360,:,:) = [];
save('F:\My research\case2\JGR\1_Figures\Figure2\Net_Xc_X_Xp.mat','Data_CNP_map')
%load 'F:\My research\case2\JGR\Working\step3_reNPPtuaE_fct\NetCNP_Xp.mat'
% global map of net change of C storage, C storage capacity and the magnitude of disequilibrium
load('F:\My research\case2\JGR\1_Figures\Figure2\map_X12.mat')
figure
set(gcf,'position',[100 20 900 600])
maps_CNP = tight_subplot(3,3,[-0.02 -0.07],[0.035 0.001],[-0.03 0.11])
sub_lat = tight_subplot(3,1,[0.08 1],[0.088 0.049],[0.88 0.01])

Labels = {'(a)','(b)','(c)','(d)',...
          '(e)','(f)','(g)','(h)',...
          '(i)','(j)','(k)','(l)',...
          '(m)','(n)','(o)','(p)',...
          '(q)','(r)','(s)','(t)'};          
for i= 1:3
    i
    
    map_i = flipud(Data_CNP_map(:,:,i));
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    axes(maps_CNP(i))
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(map_X12)
    caxis([-10 10])
    set(gca,'box','off')
    axis off
    colorbar('off')
    %text(-0.02352,1.55, Models6{i},'HorizontalAlignment','center',...
    %    'FontName','Arial','FontSize',10);
    text(-2.649,1.9, Labels{i},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

h = colorbar
h.Location = 'southoutside'
h.Position = [0.1942,0.698,0.469,0.0161];
h.Label.FontName = 'Arial'
h.Label.FontSize = 9;
text(-6.5731,-2.1448, '\DeltaXp (KgC m^-^2)','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',10);

% Zonal mean plots
axes(sub_lat(1))
hold on
C_pot_lat = nanmean(Data_CNP_map(:,:,1),2);
CN_pot_lat = nanmean(Data_CNP_map(:,:,2),2);
CNP_pot_lat = nanmean(Data_CNP_map(:,:,3),2);

lat = 90:-0.5:-59.5;
plot(C_pot_lat,lat,'LineStyle','-','LineWidth',1.7,'Color', [0.34 0.33 0.33])
plot(CN_pot_lat,lat,'LineStyle','-','LineWidth',1.7,'Color', [1 0 0])
plot(CNP_pot_lat,lat,'LineStyle','-','LineWidth',1.7,'Color', [0 0 1])


plot([6,10],[-18,-18],'LineStyle','-','LineWidth',1.7,'Color', [0.34 0.33 0.33])
text(11,-16, 'C','HorizontalAlignment','left','FontName','Arial','FontSize',8);
plot([6,10],[-28,-28],'LineStyle','-','LineWidth',1.7,'Color', [1 0 0])
text(11,-28, 'CN','HorizontalAlignment','left','FontName','Arial','FontSize',8);
plot([6,10],[-38,-38],'LineStyle','-','LineWidth',1.7,'Color', [0 0 1])
text(11,-40, 'CNP','HorizontalAlignment','left','FontName','Arial','FontSize',8);

set(gca,'lineWidth',1.2,'box', 'on')
set(gca,'XLim',[-0.2 18],'YLim',[-60 90])
set(gca,'Fontname','Arial','FontSize',8);
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xlabel('KgC m^-^2','Fontname','Arial','FontSize',10)
plot([-0.2 18],[0 0],'k--','LineWidth',1.1)
text(2.5,77, '(d)','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);

for i= 4:6
    i
    
    map_i = flipud(Data_CNP_map(:,:,i));
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    axes(maps_CNP(i))
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(map_X12)
    caxis([-20 20])
    set(gca,'box','off')
    axis off
    colorbar('off')
    %text(-0.02352,1.55, Models6{i},'HorizontalAlignment','center',...
    %    'FontName','Arial','FontSize',10);
    text(-2.649,1.9, Labels{i+1},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

h = colorbar
h.Location = 'southoutside'
h.Position = [0.1942,0.3833,0.469,0.0161];
h.Label.FontName = 'Arial'
h.Label.FontSize = 9;
text(-6.5731,-2.138, '\DeltaXc (KgC m^-^2)','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',10);

% Zonal mean plots
axes(sub_lat(2))
hold on
C_cap_lat = nanmean(Data_CNP_map(:,:,4),2);
CN_cap_lat = nanmean(Data_CNP_map(:,:,5),2);
CNP_cap_lat = nanmean(Data_CNP_map(:,:,6),2);

lat = 90:-0.5:-59.5;
plot(C_cap_lat,lat,'LineStyle','-','LineWidth',1.7,'Color', [0.34 0.33 0.33])
plot(CN_cap_lat,lat,'LineStyle','-','LineWidth',1.7,'Color', [1 0 0])
plot(CNP_cap_lat,lat,'LineStyle','-','LineWidth',1.7,'Color', [0 0 1])

set(gca,'lineWidth',1.2,'box', 'on')
set(gca,'XLim',[-0.2 22],'YLim',[-60 90])
set(gca,'Fontname','Arial','FontSize',8);
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');

yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xlabel('KgC m^-^2','Fontname','Arial','FontSize',10)
plot([-0.2 22],[0 0],'k--','LineWidth',1.1)
text(3.5,77, '(h)','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);


for i= 7:9
    i
    
    map_i = flipud(Data_CNP_map(:,:,i));
    raster_map = georasterref('RasterSize',size(map_i),'Latlim',[-60 90],'Lonlim',[-180 180]);
    
    axes(maps_CNP(i))
    hold on
    axesm miller
    setm(gca,'MapLatLimit',[-60 90])
    framem('FLineWidth',1)
    framem('off')
    geoshow('landareas.shp','FaceColor','none')
    framem('FLineWidth',1)
    geoshow(map_i,raster_map, 'DisplayType','surface','Zdata',zeros(size(map_i)),'CData',map_i);
    colormap(map_X12)
    caxis([-10 10])
    set(gca,'box','off')
    axis off
    colorbar('off')
    %text(-0.02352,1.55, Models6{i},'HorizontalAlignment','center',...
    %    'FontName','Arial','FontSize',10);
    text(-2.649,1.9, Labels{i+2},'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);
    
end

h = colorbar
h.Location = 'southoutside'
h.Position = [0.1942,0.0666,0.469,0.0161];
h.Label.FontName = 'Arial'
h.Label.FontSize = 9;
text(-6.5731,-2.2449, '\DeltaX (KgC m^-^2)','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',10);

% Zonal mean plots
axes(sub_lat(3))
hold on
C_pol_lat = nanmean(Data_CNP_map(:,:,7),2);
CN_pol_lat = nanmean(Data_CNP_map(:,:,8),2);
CNP_pol_lat = nanmean(Data_CNP_map(:,:,9),2);

plot(C_pol_lat,lat,'LineStyle','-','LineWidth',1.7,'Color', [0.34 0.33 0.33])
plot(CN_pol_lat,lat,'LineStyle','-','LineWidth',1.7,'Color', [1 0 0])
plot(CNP_pol_lat,lat,'LineStyle','-','LineWidth',1.7,'Color', [0 0 1])

set(gca,'lineWidth',1.2,'box', 'on')
set(gca,'XLim',[-0.2 8.5],'YLim',[-60 90])
set(gca,'Fontname','Arial','FontSize',8);
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');

yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xlabel('KgC m^-^2','Fontname','Arial','FontSize',10)
plot([-0.2 8.5],[0 0],'k--','LineWidth',1.1)
text(1.2,77, '(l)','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);



  


    


