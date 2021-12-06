% Global estimation of 
% re_NPP * tuaE0, re_tuaE * NPP0,  re_NPP * re_tuaE and Xp
clear
clc
load('G:\My research\case2\JGR\Working\step4_PlantSoil\CNPps_NPPtuaE_maps.mat');

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

%% Figure for re_NPP * tuaE0, re_tuaE * NPP0,  re_NPP * re_tuaE and Xp 
% for plant
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


Datap_reNPPtuaE_PgG = [];
Datas_reNPPtuaE_PgG = [];

for i = 1:12
    i
    reNPPtuaE_gbP = Data_reNPPtuaE_p(:,:,i).*cellarea .* 10^6./10^15;     % convert unit from gC m-2 into Pg C 
    reNPPtuaE_gbS = Data_reNPPtuaE_s(:,:,i).*cellarea .* 10^6./10^15;
    
    Datap_reNPPtuaE_PgG(i) = nansum(reNPPtuaE_gbP(:));
    Datas_reNPPtuaE_PgG(i) = nansum(reNPPtuaE_gbS(:));
end

% plant
C3_reNPPtuaE_p = Datap_reNPPtuaE_PgG(1:3);
CN3_reNPPtuaE_p = Datap_reNPPtuaE_PgG(5:7);
CNP3_reNPPtuaE_p = Datap_reNPPtuaE_PgG(9:11);

BarData_p(1,1:3) = C3_reNPPtuaE_p;
BarData_p(2,1:3) = CN3_reNPPtuaE_p;
BarData_p(3,1:3) = CNP3_reNPPtuaE_p;

Pt_Xc = nansum(BarData_p,2);
pt_Xp = [Datap_reNPPtuaE_PgG(4), Datap_reNPPtuaE_PgG(8), Datap_reNPPtuaE_PgG(12)]';
Pt_X = Pt_Xc - pt_Xp;
BarData_p_XXp = [Pt_X pt_Xp];

% soil
C3_reNPPtuaE_s = Datas_reNPPtuaE_PgG(1:3);
CN3_reNPPtuaE_s = Datas_reNPPtuaE_PgG(5:7);
CNP3_reNPPtuaE_s = Datas_reNPPtuaE_PgG(9:11);

BarData_s(1,1:3) = C3_reNPPtuaE_s;
BarData_s(2,1:3) = CN3_reNPPtuaE_s;
BarData_s(3,1:3) = CNP3_reNPPtuaE_s;

Sl_Xc = nansum(BarData_s,2);
Sl_Xp = [Datas_reNPPtuaE_PgG(4), Datas_reNPPtuaE_PgG(8), Datas_reNPPtuaE_PgG(12)]';
Sl_X = Sl_Xc - Sl_Xp;

BarData_s_XXp = [Sl_X Sl_Xp];



% plant
figure
set(gcf,'position',[100 100 970 420])
CNP_fig = tight_subplot(1,2,[0 0.08],[0.15 0.02],[0.1 0.02])

axes(CNP_fig(1));
hold on
Xplant=[1,4,7]
pt_bar = bar(Xplant, BarData_p,'stacked','BarWidth',0.2)
pt_bar(1).FaceColor = [0.47,0.67,0.19];
pt_bar(2).FaceColor = [0.93,0.69,0.13];
pt_bar(3).FaceColor = [1.00,0.45,0.00];
pt_bar(1).EdgeColor  = 'none';
pt_bar(2).EdgeColor  = 'none';
pt_bar(3).EdgeColor  = 'none';

Xplant2=[1.7,4.7,7.7]
pt_XXp_bar = bar(Xplant2, BarData_p_XXp,'stacked','BarWidth',0.2)
pt_XXp_bar(1).FaceColor = [0.25,0.25,0.99];
pt_XXp_bar(2).FaceColor = "none";
pt_XXp_bar(1).EdgeColor = [0.00,0.00,1.00];
pt_XXp_bar(1).LineWidth = 1.3;
pt_XXp_bar(2).EdgeColor = [0.00,0.00,1.00];
pt_XXp_bar(2).LineWidth = 1.3;


set(gca,'linewidth',1.2,'box','on')
set(gca,'XLim',[-0.4 8.5],'YLim',[0 200]);
set(gca,'Fontname','Arial','FontSize',12);
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');

legCNP = legend({'\tau_0 X \DeltaCinput','\Delta\tau X Cinput_0 ','\Delta\tau X \DeltaCinput','\DeltaX','\DeltaXp'});
set(legCNP, 'interpreter', 'tex')
set(legCNP,'color','none','EdgeColor','none','Fontname','Arial','Fontsize',11)

ylabel('Plant C change (PgC)','Fontname','Arial','FontSize',14)
xticks([1.3, 4.3, 7.3]);
xticklabels({'C-only','CN','CNP'})
text(-0.1, 188,'(a)','Fontname','Arial','FontSize',14)


% soil
axes(CNP_fig(2));
hold on

Xsoil=[1,4,7]
sl_bar = bar(Xsoil, BarData_s,'stacked','BarWidth',0.2)
sl_bar(1).FaceColor = [0.47,0.67,0.19];
sl_bar(2).FaceColor = [0.93,0.69,0.13];
sl_bar(3).FaceColor = [1.00,0.45,0.00];
sl_bar(1).EdgeColor  = 'none';
sl_bar(2).EdgeColor  = 'none';
sl_bar(3).EdgeColor  = 'none';

Xsoil_XXp = [1.7,4.7,7.7];
sl_XXp_bar = bar(Xsoil_XXp, BarData_s_XXp,'stacked','BarWidth',0.2)
sl_XXp_bar(1).FaceColor = [0.25,0.25,0.99];
sl_XXp_bar(2).FaceColor = "none";
sl_XXp_bar(1).EdgeColor = [0.00,0.00,1.00];
sl_XXp_bar(1).LineWidth = 1.3;
sl_XXp_bar(2).EdgeColor = [0.00,0.00,1.00];
sl_XXp_bar(2).LineWidth = 1.3;

set(gca,'linewidth',1.2,'box','on')
set(gca,'XLim',[-0.4 8.5],'YLim',[0 580]);
set(gca,'Fontname','Arial','FontSize',12);
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');

ylabel('Soil C change (PgC)','Fontname','Arial','FontSize',14)
xticks([1.3, 4.3, 7.3]);
xticklabels({'C-only','CN','CNP'})
text(-0.1, 545,'(b)','Fontname','Arial','FontSize',14)


CNP_sl_fig = tight_subplot(1,1,[0 0.08],[0.65 0.06],[0.76 0.04])
axes(CNP_sl_fig);
hold on

C3_reNPPtuaE23_s = Datas_reNPPtuaE_PgG(2:3);
CN3_reNPPtuaE23_s = Datas_reNPPtuaE_PgG(6:7);
CNP3_reNPPtuaE23_s = Datas_reNPPtuaE_PgG(10:11);

BarData23_s(1,1:2) = C3_reNPPtuaE23_s;
BarData23_s(2,1:2) = CN3_reNPPtuaE23_s;
BarData23_s(3,1:2) = CNP3_reNPPtuaE23_s;

Xsoil=[1,2,3]
sl_bar = bar(Xsoil, BarData23_s,'stacked','BarWidth',0.5)
sl_bar(1).FaceColor = [0.93,0.69,0.13];
sl_bar(2).FaceColor = [1.00,0.45,0.00];
sl_bar(1).EdgeColor  = 'none';
sl_bar(2).EdgeColor  = 'none';


set(gca,'linewidth',1.2,'box','on')
set(gca,'XLim',[0.3 3.7],'YLim',[-1.5 3]);
set(gca,'Fontname','Arial','FontSize',12);

xticks([1, 2, 3]);
xticklabels({'C-only','CN','CNP'})
ylabel('Soil C change (PgC)','Fontname','Arial','FontSize',12)
text(3.1, 2.4,'(c)','Fontname','Arial','FontSize',14)

