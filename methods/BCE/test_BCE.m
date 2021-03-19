% Test script for DICLENS family

scoreExt = {'NMI'};

%% Load data
[data,target] = pplk_loadData('iris');
%[data,target] = pplk_loadData('real\iris_orig');
[N D] = size(data);
%% Generate ensemble
params = pplk_setParamsDefault();
%params.KM_nRuns = 10;
ensSize = 10;
K = floor(sqrt(N));

labelsEns = pplk_genEns(data,{'KM',ensSize,K,'rand'},params);

Kcons = max(target);


% params

Ki = zeros(1,ensSize);
for iL = 1:ensSize
    lbl = labelsEns(:,iL);
    Ki(iL) = length(unique(lbl));
end
alphaInit = rand(Kcons,1);
betaInit = cell(1,ensSize);
for i=1:length(betaInit)
    temp = rand(Kcons,Ki(i));
    temp = temp./(sum(temp,2)*ones(1,Ki(i)));
    betaInit{i} = temp;
end

params.BCE_alphaInit = alphaInit;
params.BCE_betaInit = betaInit;

%% fresh & faster implementation
ticID2 = tic();
[labelsCons1,numClust]= pplk_consEns(labelsEns,Kcons,'BCE',params);
t2 = toc(ticID2);
%pplk_scatterPlot(data,labelsCons1);
[~,score2] = pplk_validExt(target,labelsCons1,scoreExt);
fprintf(1,'BCE-fast | score: %f, time: %f\n', score2, t2);

%% (original implementation)
ticID1 = tic();
[labelsCons2,numClust] = pplk_consEns(labelsEns,Kcons,'BCE_ORIG',params);
t1 = toc(ticID1);

%pplk_scatterPlot(data,labelsCons2);
[~,score1] = pplk_validExt(target,labelsCons2,scoreExt);
fprintf(1,'BCE-ORIG | score: %f, time: %f\n', score1, t1);

isequal(labelsCons1, labelsCons2)
