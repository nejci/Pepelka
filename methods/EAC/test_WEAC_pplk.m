% test script for WEAC consensus algorithm
% reproducing results from table X in Duarte et al., 2005 (script was
% modified since)

eCVI = {'ARI'};
%% Load data
%[data,target] = pplk_loadData('artificial\spiral');
%[data,target] = pplk_loadData('artificial\halfring');
%[data,target] = pplk_loadData('artificial\ring');
%[data,target] = pplk_loadData('real\iris_orig');
%[data,target] = pplk_loadData('real\wine');
[data,target] = pplk_loadData('real\iris'); % try also iris_orig


%% Generate ensemble
params = pplk_setParamsDefault();

%labelsEns = pplk_genEns(data,{'SL',1,[10,30],'fixed'}, params);
%labelsEns = pplk_genEns(data,{'AL',1,[10,30],'fixed'}, params);
%labelsEns = pplk_genEns(data,{'CL',1,[10,30],'fixed'}, params);
labelsEns = pplk_genEns(data,{'KM',200,[10,30],'rand'},params);

Kcons = max(target);

%% consensus using EAC-SL
[labelsCons]= pplk_consEns(labelsEns,Kcons,'EAC-SL');

Kout_EAC_SL = length(unique(labelsCons));
score_EAC_SL = pplk_validExt(target,labelsCons,eCVI)
pplk_scatterPlot(data,labelsCons);
%% consensus using EAC-AL
[labelsCons]= pplk_consEns(labelsEns,Kcons,'EAC-AL');

Kout_EAC_AL = length(unique(labelsCons));
pplk_scatterPlot(data,labelsCons);
score_EAC_AL = pplk_validExt(target,labelsCons,eCVI)

%% consensus using EAC-WL
[labelsCons]= pplk_consEns(labelsEns,Kcons,'EAC-WL');

Kout_EAC_WL = length(unique(labelsCons));
pplk_scatterPlot(data,labelsCons);
score_EAC_WL = pplk_validExt(target,labelsCons,eCVI)

%% WEAC with Pepelka wrapper

params.WEAC_CVI = {'DN','DB','SDBW','CH','SIL','I','XB'};
params.WEAC_data = data;
params.WEAC_CVImat = [];
params.WEAC_unifyMeth = 'minmax';
params.WEAC_reduceMeth = 'SPEC';
params.WEAC_reduceDim = 'DANCo';
params.WEAC_weightMeth = 'wRankAggreg2';
params.WEAC_weightMode = 'stuart';

weights = ones(numEns,1)*0;
params.WEAC_weights = weights;

[labelsCons,numClust]= pplk_consEns(labelsEns,Kcons,'WEAC-SL',params);
score_WEAC = pplk_validExt(target,labelsCons,eCVI)

% consensus using WEAC
%% compute CVI
CVI = {'DN','DB','SDBW','CH','SIL','I','XB'};
%CVI = {'DN','DB','SDBW','CH','SIL','I','XB','XB*','DB*','DNg','DNs','CON','SD','SF','CI','GAMMA'};

[validStrct, list] = pplk_validInt(data, labelsEns, CVI, []);
PRM_raw = list';

%% compute weights and WEAC
consClusterer = 'single';
unifyMethod = 'minmax';
reduceMethod = 'SPEC';
reduceDim = 'DANCo';
weightMethod = 'wRankAggreg2';
weightMode = 'stuart';

PRMopt = [];
PRMopt.CVImat = PRM_raw;
PRMopt.reduceDim = reduceDim;

[weights,PRM,featInd] = ...
    pplk_partitionRelevance(data,labelsEns,CVI,...
    unifyMethod,reduceMethod,weightMethod,weightMode,PRMopt);

numEns = size(labelsEns,2);
weights = ones(numEns,1)*0;

WEACopt = [];
WEACopt.weights = weights;

[labelsCons]= WEAC(data,labelsEns,Kcons,consClusterer,CVI,weightMethod,WEACopt);

score_WEAC = pplk_validExt(target,labelsCons,eCVI)


