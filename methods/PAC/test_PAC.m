% test script for PA consensus algorithm


%% Load data
%[data,target] = pplk_loadData('artificial\halfring');
[data,target] = pplk_loadData('real\iris_orig'); % try also iris_orig

%% Generate ensemble
params = pplk_setParamsDefault();
params.KM_nRuns = 1;
params.subsampling = {'cols','rand'};
ensSize = 50;
K = [2,8];

[labelsEns,moreInfo] = pplk_genEns(data,{'KM',ensSize,K,'rand'},params);

Kcons = 3;

%% consensus using PA
[labelsCons, numClust] = PAC(labelsEns,size(data,2),[],'SL');

Kout = length(unique(labelsCons))
pplk_scatterPlot(data,labelsCons);
[~,scorePA] = pplk_validExt(target,labelsCons,{'CA','NMI'})

%% consensus using EAC
[labelsCons]= pplk_consEns(labelsEns,[],'EAC-SL');

%Kout = length(unique(labelsCons))
%pplk_scatterPlot(data,labelsCons);
[~,scoreEAC] = pplk_validExt(target,labelsCons,{'CA','NMI'})

%% Strehl
[labelsCons] = pplk_consEns(labelsEns,Kout,'STREHL-MCLA');
%pplk_scatterPlot(data,labelsCons);
[~,scoreMCLA] = pplk_validExt(target,labelsCons,{'CA','NMI'})

%% HBGF
[labelsCons] = pplk_consEns(labelsEns,Kout,'HBGF');
%pplk_scatterPlot(data,labelsCons);
[~,scoreHGBF] = pplk_validExt(target,labelsCons,{'CA','NMI'})
