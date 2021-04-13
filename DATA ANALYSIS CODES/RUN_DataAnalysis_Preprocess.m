clear all
%%%
codedir = 'Z:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\CycIF Codes\Utility Functions\';
addpath(codedir)
addpath('Z:\sorger\data\RareCyte\Ana\2021-03-PDAC-Pin1\ANALYSIS\Codes Used\Pipeline')
addpath('Z:\sorger\data\RareCyte\Ana\2021-03-PDAC-Pin1\ANALYSIS\Codes Used\Pipeline\norm\')

%%%
options.date = '20210318';
filename.basefolder = 'Z:\sorger\data\RareCyte\Ana\2021-03-PDAC-Pin1\';
filename.suffix = ['_Results_' options.date '.mat'];
% filename.analfolder = [filename.basefolder 'Data Analysis\'];
filename.analfolder = [filename.basefolder 'ANALYSIS\'];
filename.resufolder = 'Analysis_Results\'; 
filename.roifolder  = 'ROIs\';
filename.montfolder = 'TumorMasks\MontageforROI_Lv3\';
filename.montsuffix = '_montage.tif';

filename.folders = {'PDAC_1'   ,'PDAC_2'   ,'PDAC_3'   ,'PDAC_4'   , 'PDAC_5_1', ...
                    'PDAC_5_2' ,'PDAC_6'   ,'PDAC_7'   ,'PDAC_9'   , 'PDAC_10'};   
filename.tissues = {'PDAC_1'   ,'PDAC_2'   ,'PDAC_3'   ,'PDAC_4'   , 'PDAC_5_1'  ,...
                    'PDAC_5_2' ,'PDAC_6'   ,'PDAC_7'   ,'PDAC_9'   , 'PDAC_10'};   


options.Markers =  { 	
    'DAPI0', 'BKG488'   , 'BKG555'   , 'BKG647' , ...
    'DAPI1', 'BKG488'   , 'BKG555'   , 'BKG647' , ...
    'DAPI2', 'BKG488'   , 'CD73'     , 'PIN1'   , ...
    'DAPI3', 'DPB1'     , 'PDGFRa'   , 'FAP'    , ...
    'DAPI3', 'Vimentin' , 'CD90'     , 'PIN1-sc', ...
    'DAPI3', 'panCK'    , 'CD44'     , 'CD45'   , ...
    'DAPI3', 'KI67'     , 'aSMA'     , 'BKG647' , ...
    'DAPI3', 'CD4'      , 'CD3d'     , 'CD8a'   , ...
    'DAPI3', 'CD163'    , 'CD68'     , 'CD20'   , ...
};
options.maxround = length(options.Markers);
options.magnification = 20;
options.FigOpt = 0;
% options.date = '20201009';

mkdir([filename.analfolder filename.resufolder])
save([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'],'filename','options')

%% Step 1 aggregate tissues together
clearvars -except filename options
[AggrResults, MorpResults] = PreProcess_Step1_Aggregation_v2(filename, options, options.FigOpt);

% save
save([filename.analfolder filename.resufolder 'Results_Aggr_' options.date '.mat'],'AggrResults');
save([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'],'MorpResults');
save([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'],'options', 'filename');

%% Step 2 - Filter the data
clearvars -except filename options
load([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Aggr_' options.date '.mat'])
options.FigOpt = 0;
%options.maxround = 36;

% STEP 2: additional parameters for filtering          
options.Filtering.folder = [filename.analfolder filename.resufolder 'Step2_'];
options.Filtering.Index_Names = filename.tissues;
options.Filtering.thresholds.foldDAPI_th = 0.6; %0.5 
options.Filtering.thresholds.absDAPI_th = 8; % 9 
options.Filtering.thresholds.solidity = 0.8; 
options.Filtering.thresholds.area_low = 20;    
options.Filtering.thresholds.area_high = 700; 
options.Filtering.maxround = options.maxround;

%%
[Filter,report] = PreProcess_Step2_Filter(AggrResults, options.Filtering, options.FigOpt);
save([filename.analfolder filename.resufolder 'Results_Filt_' options.date '.mat'],'Filter','report')

%% Step 3 - Normalize data
% close all
clearvars -except filename options
load([filename.analfolder filename.resufolder 'Results_Filt_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Aggr_' options.date '.mat'])
load([filename.analfolder filename.resufolder 'Results_Morp_' options.date '.mat'])

rng(11)
close all

%%
% STEP 3: additional parameters for normalization
options.Norm.Reps = 5;
options.Norm.FigSettings.FigFlag = 1; % to save
options.Norm.FigSettings.Folder = [filename.analfolder filename.resufolder 'Step3_NormPrints\' ];
options.Norm.FigSettings.Debug = 0; % to view
options.Norm.Channels = setdiff(1:length(options.Markers), 1:4:length(options.Markers));  % non-dapi channels
% options.Norm.Channels = [14];
% default to zeros
options.Norm.Priors = zeros(length(options.Markers),1); 
options.Norm.OverExpr = zeros(length(options.Markers),1);
options.Norm.FinalEst = zeros(length(options.Markers),1);


priors = zeros(length(options.Markers),1);
priors([11 12 19 20 27 30 31 34]) = +1;
options.Norm.Priors = priors;

overexpr = [0    0     0    0      ...  'DAPI0',  'BKG 488' ,   'BKG 555',   'BKG 647', ...
            0    0     0    0      ...  'DAPI1',  'BKG 488' ,   'BKG 555',   'BKG 647', ...
            0    0     0.08 0.66   ...  'DAPI2',  'BKG 488' ,   'CD73'   ,   'PIN1'   , ...
            0    0.15     0    0   ...  'DAPI3',  'MHC-II'  ,   'PDGFRa' ,   'FAP'    , ...
            0    0     0.1  0      ...  'DAPI4',  'Vimentin',   'CD90'   ,   'PIN1sc' , ...
            0    0     0    0      ...  'DAPI5',  'panCK'   ,   'CD44'   ,   'CD45'   , ...
            0    0     0.08 0      ...  'DAPI6',  'Ki67'    ,   'aSMA'   ,   'BKG647' , ...
            0    0     0    0      ...  'DAPI7',  'CD4'     ,   'CD3d'   ,   'CD8a'   , ...
            0    0     0    0      ...  'DAPI8',  'CD163'   ,   'CD68'   ,   'CD20'   , ...
           ]; 

options.Norm.OverExpr = overexpr;

filter = Filter.all;
data_nuc = log2(double(AggrResults.MedianNucSign)+1);
data_cyt = log2(double(AggrResults.MedianCytSign)+1);
options.Norm.CellNum = min([50000, length(AggrResults.MedianNucSign)]);

% normalize
options.Norm.IsNuc = 1;
[NucNorm, cutoffs_nuc, mults_nuc] = norm_main(data_nuc, filter, options.Markers, options.Norm);
options.Norm.IsNuc = 0;
[CytNorm, cutoffs_cyt, mults_cyt] = norm_main(data_cyt, filter, options.Markers, options.Norm);

% save
redmap = [linspace(0,255,128) zeros(1,128)+255 ]/255;
blumap = [zeros(1,128)+255 flip(linspace(0,255,128))]/255;
gremap = [linspace(128,255,128) flip(linspace(128,255,128))]/255;
NormResults.colorMap = [redmap' gremap' blumap'];
NormResults.CellID = (1:size(AggrResults.MeanNucSign,1))'; 

NormResults.MedianNucNorm = int16(round(NucNorm*1000,0));
NormResults.MedianCytNorm = int16(round(CytNorm*1000,0));
NormResults.nuc_add_fact = int16(1000.*cutoffs_nuc); 
NormResults.nuc_mult_fact = int16(1000.*mults_nuc);
NormResults.cyt_add_fact = int16(1000.*cutoffs_cyt); 
NormResults.cyt_mult_fact = int16(1000.*mults_cyt);


save([filename.analfolder filename.resufolder 'Results_Norm_' options.date '.mat'],'NormResults');
save([filename.analfolder filename.resufolder 'Results_Settings_' options.date '.mat'],'filename','options');

%% 
idcs = find(Filter.all(:,options.maxround) == 1);
rand_idcs = randsample(length(idcs),100000);
figure
dapi = 1:4:33;
markers = ~ismember(10:36,dapi);
vec = 10:36;
markers = vec(markers);
counter = 1;
myaxes = [];
for i = markers
    data_nuc = NormResults.MedianNucNorm(idcs(rand_idcs),i);
    data_cyt = NormResults.MedianCytNorm(idcs(rand_idcs),i);
    subplot(3,7,counter)
    myaxes(counter) = gca;
    ksdensity(data_nuc) %, 'Bandwidth', 0.001);
    hold on
    ksdensity(data_cyt) %, 'Bandwidth', 0.001)
    hold on
    xline(0, '--r');
    title([num2str(i) ' - ' options.Markers{i}])
%     title(options.Markers{i})

    disp(options.Markers{i})
    counter = counter + 1;
end

for i = 1:length(myaxes)
   xlim(myaxes(i),[-3000,3000])
   ylim(myaxes(i),[0,0.0018])
   set(myaxes(i),'FontSize',20)
end
%%
ch = 27;
idcs = find(Filter.all(:,options.maxround) == 1);
rand_idcs = randsample(length(idcs),100000);

figure
data_nuc = NormResults.MedianNucNorm(idcs(rand_idcs),ch);
data_cyt = NormResults.MedianCytNorm(idcs(rand_idcs),ch);
ksdensity(data_nuc) %, 'Bandwidth', 0.001);
hold on
ksdensity(data_cyt) %, 'Bandwidth', 0.001)
hold on
xline(0, '--r');
title([num2str(i) ' - ' options.Markers{ch}])
disp(options.Markers{ch})