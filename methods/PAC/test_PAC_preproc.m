% test script for PA consensus algorithm
% using pre/post processing
% I dont see any improvement, it is actually worse.


%% Load data
%[data,target] = pplk_loadData('artificial\halfring');
[data,target] = pplk_loadData('real\iris_orig'); % try also iris_orig

% figure(); pplk_scatterPlot(data);

% normalize data
dataN = pplk_normalize(data,'zscore');

%figure(); pplk_scatterPlot(dataN);


% reflect data on a grid
epsilon = 0.1;
dataR = round(dataN./epsilon).* epsilon;
% remove duplicates
[dataR, ia, ib] = unique(dataR,'rows');

% figure(); pplk_scatterPlot(dataR,[]);

%% Generate ensemble
params = pplk_setParamsDefault();
params.KM_nRuns = 1;
ensSize = 10;
K = [10,30];

labelsEns = pplk_genEns(dataR,{'KM',ensSize,K,'rand'},params);

Kcons = 3;

%% consensus using PA
[labelsConsR, numClust] = PAC(labelsEns,size(data,2),[],'SL');

% reflect labels back to original data
labelsCons = labelsConsR(ib);

Kout = length(unique(labelsCons))
%pplk_scatterPlot(data,labelsCons);
[~,scorePA] = pplk_validExt(target,labelsCons,{'CA','NMI'})

%% consensus using EAC
[labelsConsR]= pplk_consEns(labelsEns,[],'EAC-SL');

% reflect labels back to original data
labelsCons = labelsConsR(ib);

%Kout = length(unique(labelsCons))
%pplk_scatterPlot(data,labelsCons);
[~,scoreEAC] = pplk_validExt(target,labelsCons,{'CA','NMI'})

%% Strehl
[labelsConsR] = pplk_consEns(labelsEns,Kout,'STREHL-MCLA');
% reflect labels back to original data
labelsCons = labelsConsR(ib);
%pplk_scatterPlot(data,labelsCons);
[~,scoreMCLA] = pplk_validExt(target,labelsCons,{'CA','NMI'})

%% HBGF
[labelsConsR] = pplk_consEns(labelsEns,Kout,'HBGF');
% reflect labels back to original data
labelsCons = labelsConsR(ib);
%pplk_scatterPlot(data,labelsCons);
[~,scoreHGBF] = pplk_validExt(target,labelsCons,{'CA','NMI'})
