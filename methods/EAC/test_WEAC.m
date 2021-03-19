% test script for WEAC consensus algorithm
% reproducing results from table X in Duarte et al., 2005

%% Load data
[data,target] = pplk_loadData('artificial\spiral');
%[data,target] = pplk_loadData('artificial\halfring');
%[data,target] = pplk_loadData('artificial\ring');
%[data,target] = pplk_loadData('real\iris_orig');
%[data,target] = pplk_loadData('real\wine');
%[data,target] = pplk_loadData('real\iris'); % try also iris_orig


%% Generate ensemble
params = pplk_setParamsDefault();

labelsEns = pplk_genEns(data,{'SL',1,[10,30],'fixed'}, params);
%labelsEns = pplk_genEns(data,{'AL',1,[10,30],'fixed'}, params);
%labelsEns = pplk_genEns(data,{'CL',1,[10,30],'fixed'}, params);
%labelsEns = pplk_genEns(data,{'KM',200,[10,30],'rand'},params);

Kcons = max(target);

%% consensus using EAC-SL
[labelsCons]= pplk_consEns(labelsEns,Kcons,'EAC-SL');

Kout_EAC_SL = length(unique(labelsCons));
pplk_scatterPlot(data,labelsCons);
score_EAC_SL = pplk_validExt(target,labelsCons,{'CA','NMI'})

%% consensus using EAC-AL
[labelsCons]= pplk_consEns(labelsEns,Kcons,'EAC-AL');

Kout_EAC_AL = length(unique(labelsCons));
pplk_scatterPlot(data,labelsCons);
score_EAC_AL = pplk_validExt(target,labelsCons,{'CA','NMI'})

%% consensus using EAC-WL
[labelsCons]= pplk_consEns(labelsEns,Kcons,'EAC-WL');

Kout_EAC_WL = length(unique(labelsCons));
pplk_scatterPlot(data,labelsCons);
score_EAC_WL = pplk_validExt(target,labelsCons,{'CA','NMI'})

% consensus using WEAC
%% compute CVI
CVI = {'DN','DB','SDBW','CH','SIL','I','XB'};
%CVI = {'DN','DB','SDBW','CH','SIL','I','XB','XB*','DB*','DNg','DNs','CON','SD','SF','CI','GAMMA'};

[validStrct, list] = pplk_validInt(data, labelsEns, CVI, []);
PRM_raw = list';

%% compute weights and WEAC
consClusterer = 'SL';
unifyMethod = 'minmax';
reduceMethod = 'none'; % method for dimensionality reduction
reduceDim = 'DANCo'; % number of remaining features | name of intrinsic dimensionality estimator
weightMethod = 'wMean2';
weightMode = 'CLK';

PRMopt = [];
PRMopt.CVImat = PRM_raw;
PRMopt.reduceDim = reduceDim;

[weights,PRM,featInd] = ...
    pplk_partitionRelevance(data,labelsEns,CVI,...
    unifyMethod,reduceMethod,weightMethod,weightMode,PRMopt);

WEACopt = [];
WEACopt.weights = weights;

[labelsCons]= WEAC(data,labelsEns,Kcons,consClusterer,CVI,weightMethod,WEACopt);

numEns = size(labelsCons,2);
eCVI = {'CA'};
score_WEAC = zeros(length(eCVI),numEns);
for i = 1:numEns
    Kout = length(unique(labelsCons(:,i)));
    [~,score_WEAC(:,i)]=pplk_validExt(target,labelsCons(:,i),eCVI);
end
score_WEAC
max(score_WEAC(1,:))