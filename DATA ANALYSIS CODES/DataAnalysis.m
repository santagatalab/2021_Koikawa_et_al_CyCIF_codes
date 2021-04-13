%% Loading files
clear all

addpath('Z:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\CycIF Codes\Utility Functions');
basefolder = 'Z:\sorger\data\RareCyte\Ana\2021-03-PDAC-Pin1\';
analfolder = [basefolder 'ANALYSIS\'];
resufolder = 'Analysis_Results\';
date = '20210318';

load([ analfolder resufolder 'Results_Aggr_' date '.mat'])
load([ analfolder resufolder 'Results_Morp_' date '.mat'])
load([ analfolder resufolder 'Results_Norm_' date '.mat'])
load([ analfolder resufolder 'Results_Filt_' date '.mat'])
load([ analfolder resufolder 'Results_Settings_' date '.mat'])


filename.basefolder = basefolder;
filename.analfolder = analfolder;
filename.resufolder = resufolder;
options.date = date;
save([ analfolder resufolder 'Results_Settings_' date '.mat'], 'filename','options','-append')


%% Loading marker rounds and signals 
pinround = find(strcmp(options.Markers,'PIN1'));       %PIN1 from Gerburg lab
mhc2round = find(strcmp(options.Markers,'DPB1'));      
pinscround = find(strcmp(options.Markers,'PIN1-sc'));  %PIN1, Santa Cruz
panround = find(strcmp(options.Markers,'panCK'));     
cd44round = find(strcmp(options.Markers,'CD44'));      
cd45round = find(strcmp(options.Markers,'CD45'));      
asmaround = find(strcmp(options.Markers,'aSMA'));    
cd90round = find(strcmp(options.Markers,'CD90'));
cd73round = find(strcmp(options.Markers,'CD73'));
cd4round = find(strcmp(options.Markers,'CD4'));
cd8round = find(strcmp(options.Markers,'CD8a'));
cd20round = find(strcmp(options.Markers,'CD20'));

ind = prod(Filter.all(:,:),2) == 1;

%%Normalized Signals 
PIN1 = NormResults.MedianNucNorm(:,pinround);
MHCII = NormResults.MedianNucNorm(:,mhc2round);
PIN1SC = NormResults.MedianNucNorm(:,pinscround);
PANCK = NormResults.MedianNucNorm(:,panround);
CD44 = NormResults.MedianNucNorm(:,cd44round);
CD45 = NormResults.MedianNucNorm(:,cd45round);
ASMA = NormResults.MedianNucNorm(:,asmaround);
CD4 = NormResults.MedianNucNorm(:,cd4round);
CD8 = NormResults.MedianNucNorm(:,cd8round);
CD20 = NormResults.MedianNucNorm(:,cd20round);

%% kmeans clustering
rng(22)
ind = prod(Filter.all(:,[asmaround cd44round cd45round cd90round cd73round panround mhc2round]),2) == 1;
ind_all = ind & CD45 < -200 & PANCK < 0;

data = [ASMA CD44 MHCII PIN1 PIN1SC];
data = double(data(ind_all,:));

labels = options.Markers([asmaround cd44round mhc2round pinround pinscround]);
numclust = 20 ;
iters = 300;
rep = 4;

cmap = NormResults.colorMap;
ax_lim = [-500 500];
 
[hier_clus,hier_idx,hier_C,hier_D,order] = cluster_cycIFdata_extradata(data(:,1:3),numclust, labels,1,iters,rep,ax_lim,cmap,'',data(:,4:5),'');


extradata = data(:,4:5);
hier_clus_extra = zeros(max(hier_idx),size(extradata,2));
for i = 1:max(hier_idx)
        hier_clus_extra(i,:) = mean(extradata(find(hier_idx==i),:),1);
end 

clus_reordered = zeros(size(hier_clus));
clus_extra_reordered = zeros(size(hier_clus_extra));
reorder = [14 11 18 19 16 15 17 1 2 3 4 5 6 7 8 9 10 12 13 20]; %sort clusters based on asma, cd44, and dpb1 expression levels
counter = 1;

for i = reorder
   clus_reordered(counter,:) = hier_clus(i,:);
   clus_extra_reordered(counter,:) = hier_clus_extra(i,:);
   counter = counter + 1;
end

labels = options.Markers([asmaround cd44round mhc2round pinround pinscround]);

figure
imagesc([clus_reordered clus_extra_reordered])
colorbar
colormap(NormResults.colorMap)
caxis([-500 500])
xticks = linspace(1, numel(labels), numel(labels));
set(gca, 'XTick', xticks, 'XTickLabel', labels, 'FontWeight','bold','FontSize',20)

%% boxplots: immune, epithelial, CAFs
asma_clus = [14 11]; 
cd44_clus = 18;
dpb1_clus = 19; 

asma_idx_1 = cluster.hier_idx == asma_clus(1);
asma_idx_2 = cluster.hier_idx == asma_clus(2);
cd44_idx = cluster.hier_idx == cd44_clus;
dpb1_idx = cluster.hier_idx == dpb1_clus;
pan_ind = PANCK > 0 & CD45 < 0;
immu_ind = CD45 > 0 & PANCK < 0;

PIN1SC_new = data(:,5);
PIN1_new = data(:,4);
ASMA_new = data(:,1);

asma1_pin1 = PIN1SC_new(asma_idx_1);
asma2_pin1 = PIN1SC_new(asma_idx_2);
cd44_pin1 = PIN1SC_new(cd44_idx);
dpb1_pin1 = PIN1SC_new(dpb1_idx);
pan_pin1 = PIN1SC(pan_ind);
immu_pin1 = PIN1SC(immu_ind);

figure
ax5 = subplot(1,6,1);
boxplot(pan_pin1,'Symbol','')
title('panCK^+','FontSize',20)
ylabel('PIN1-sc','FontSize',20)
hold on
plot(mean(pan_pin1),'o','MarkerFaceColor','r')

ax6 = subplot(1,6,2);
boxplot(immu_pin1,'Symbol','')
title('CD45^+','FontSize',20)
hold on
plot(mean(immu_pin1),'o','MarkerFaceColor','r')

ax1 = subplot(1,6,3);
boxplot(asma1_pin1,'Symbol','')
title('aSMA^+^+','FontSize',20)
hold on
plot(mean(asma1_pin1),'o','MarkerFaceColor','r')
ax2 = subplot(1,6,4);
boxplot(asma2_pin1,'Symbol','')
title('aSMA^+','FontSize',20)
% yticks([])
hold on
plot(mean(asma2_pin1),'o','MarkerFaceColor','r')
ax3 = subplot(1,6,5);
boxplot(cd44_pin1,'Symbol','')
title('CD44^+','FontSize',20)
% yticks([])
hold on
plot(mean(cd44_pin1),'o','MarkerFaceColor','r')
ax4 = subplot(1,6,6);
boxplot(dpb1_pin1,'Symbol','')
title('HLA-DPB1^+','FontSize',20)
% yticks([])
hold on
plot(mean(dpb1_pin1),'o','MarkerFaceColor','r')
linkaxes([ax1 ax2 ax3 ax4 ax5 ax6])
ylim([-7000 4000])

%% boxplots: T cells 
cd8_ind = PANCK < 0 & CD45 > 0 & CD8 > 0 & CD4 < 0;
cd4_ind = PANCK < 0 & CD45 > 0 & CD4 > 0 & CD8 < 0;

cd8_pin1 = PIN1SC(cd8_ind);
cd4_pin1 = PIN1SC(cd4_ind);

figure 
ax1 = subplot(1,3,1);
boxplot(immu_pin1,'Symbol','')
title('CD45+','FontSize',20)
ylabel('PIN1-sc','FontSize',20)
hold on
plot(mean(immu_pin1),'o','MarkerFaceColor','r')

ax2 = subplot(1,3,2);
boxplot(cd8_pin1,'Symbol','')
title('CD8^+CD4^-','FontSize',20)
hold on
plot(mean(cd8_pin1),'o','MarkerFaceColor','r')

ax3 = subplot(1,3,3);
boxplot(cd4_pin1,'Symbol','')
title('CD4^+CD8^-','FontSize',20)
hold on
plot(mean(cd4_pin1),'o','MarkerFaceColor','r')

linkaxes([ax1 ax2 ax3])
ylim([-2500 3000])
