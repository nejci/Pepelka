% test script for EAC consensus algorithm
% reproducing results from table 3 in Fred & Jain, 2005

%% Load data
%[data,target] = pplk_loadData('artificial\halfring');
[data,target] = pplk_loadData('real\iris'); % try also iris_orig

%% Generate ensemble
params = pplk_setParamsDefault();
params.KM_nRuns = 10;
ensSize = 50;
K = 3;

labelsEns = pplk_genEns(data,{'KM',ensSize,K,'fixed'},params);

Kcons = 3;

%% consensus using EAC
[labelsCons]= pplk_consEns(labelsEns,Kcons,'EAC-SL');

Kout = length(unique(labelsCons))
pplk_scatterPlot(data,labelsCons);
pplk_validExt(target,labelsCons,{'CA','NMI'})

%% consensus using Strehl's method
[labelsCons2, Kout2] = pplk_consEns(labelsEns,Kcons,'STREHL-MCLA');
pplk_scatterPlot(data,labelsCons2);
pplk_validExt(target,labelsCons2,{'CA','NMI'})