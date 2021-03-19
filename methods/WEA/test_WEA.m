% test script for WEA consensus algorithm


%% Load data
%[data,target] = pplk_loadData('artificial\halfring');
[data,target] = pplk_loadData('real\iris_orig'); % try also iris_orig

%% Generate ensemble
params = pplk_setParamsDefault();
params.KM_nRuns = 1;
params.subsampling = {'none'}; %{'cols','rand'};
%params.subsampling = {'cols','rand'};
ensSize = 50;
K = [2,8];

[labelsEns,moreInfo] = pplk_genEns(data,{'KM',ensSize,K,'rand'},params);

Kcons = 3;

%% consensus using WEA
%[labelsCons, numClust] = WEA(labelsEns,data,'data','euclidean',[],'SL',1);
paramsWEA = [];
paramsWEA.WEA_data = data;
paramsWEA.WEA_dataMode = 'data';
paramsWEA.WEA_dataDist = 'euclidean';
paramsWEA.WEA_normalize = 1;

[labelsCons]= pplk_consEns(labelsEns,[],'WEA-SL', paramsWEA);

Kout = length(unique(labelsCons))
pplk_scatterPlot(data,labelsCons);
[~,scoreWEA] = pplk_validExt(target,labelsCons,{'CA','NMI'})

%% consensus using PAC
paramsPAC = [];
paramsPAC.PAC_dim = size(data,2);
[labelsCons]= pplk_consEns(labelsEns,[],'PAC-SL', paramsPAC);

%Kout = length(unique(labelsCons))
%pplk_scatterPlot(data,labelsCons);
[~,scorePAC] = pplk_validExt(target,labelsCons,{'CA','NMI'})

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
