function [labelsCons, numClust] = WEA(labelsEns,data,dataMode,dataDist,K,consFun,normalize)
% [labelsCons, numClust] = WEA(labelsEns,dimdata,dataMode,dataDist,K,consFun)
% Weighted evidence accumulation clustering (Vega-Pons et al., 2009 and 2011)
%--------------------------------------------------------------------------
% INPUTS
% labelsEns     (matrix)  clustering ensemble; each COLUMN corresponds
%						  to one clustering (standard in Pepelka).
%
% data          (matrix)  original data from which labelsEns were generated
% dataMode      (string)  data is:
%                         'data' - data matrix [N X dim]
%                         'dist' - dissimilarity matrix [N X N]
%                         'sim'  - similarity matrix [N X N]
%
% dataDist      (scalar)  if dataMode is 'data', dataDist is used to
%                         calculate pair-wise distances between data. Can
%                         be any distance supported by pdist function.
%                         Leave empty otherwise.
%
% K             (scalar)  number of clusters in output consensus
%	                      clustering; leave empty [] to automatically
%	                      determine this number by lifetime criterion.
%
% consFun       (string)  which clustering metod to use to cluster
%                           co-occurence matrix:
%                             'single'    - nearest distance (default)
%                             'complete'  - furthest distance
%                             'average'   - unweighted average distance (UPGMA)
%                                             (also known as group average)
%                             'weighted'  - weighted average distance (WPGMA)
%                             'centroid'  - unweighted center of mass distance (UPGMC)
%                             'median'    - weighted center of mass distance (WPGMC)
%                             'ward'      - inner squared distance (min variance algorithm)
%               (cell)     combination of methods in cell string
%
% normalize     (logical)  should cluster properties be normalized? Default is 1. 
%
% OUTPUTS
% labelsCons    (vector)    data labels - consensus clustering
% numClust      (scalar)    number of clusters in consensus clustering
%--------------------------------------------------------------------------
% EXAMPLE
%     %% Load data
%     [data,target] = pplk_loadData('real\iris_orig');
% 
%     %% Generate ensemble
%     params = pplk_setParamsDefault();
%     params.KM_nRuns = 10;
%     ensSize = 50;
%     K = 3;
%     labelsEns = pplk_genEns(data,{'KM',ensSize,K,'fixed'},params);
% 
%     %% consensus using WEA
%     Kcons = []; % Change if you want to control the number of clusters in
%                 % the consensus partition.
%     params.WEA_data = data;
%     params.WEA_dataMode = 'data';
%     params.WEA_dataDist = 'euclidean';
%     params.WEA_normalize = 1;
% 
%     [labelsCons, numClust] = pplk_consEns(labelsEns,Kcons,'WEA-SL',params);
% 
%     pplk_scatterPlot(data,labelsCons);
%     pplk_validExt(target,labelsCons,{'CA','NMI'})
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013  Nejc Ilc
% Part of Pepelka package.
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 13-September-2013 by Nejc Ilc
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

if ~exist('K','var')
    K=[]; % use lifetime cirteria to determine number of clusters
end

if ~exist('dataMode','var')
    dataMode = [];
end
if ~exist('dataDist','var')
    dataDist = [];
end

if ~exist('consFun','var') || isempty(consFun)
    consFun = 'single';
end

if ~exist('normalize','var')
    normalize = [];
end

% Compute co-occurence matrix, return vector of pairs
WA = computeWA(labelsEns,data,dataMode,dataDist,'vec',normalize,[]);

% Run linkage algorithm on WA vector/matrix to obtain consensus parititon.
% If K is empty, number ob final clusters is automatically determined.
params = pplk_setParamsDefault();
params.HCL_distance = [];

numConsFun = 1;
if iscell(consFun)
    numConsFun = length(consFun);
else
    consFun = {consFun};
end
labelsCons = zeros(size(labelsEns,1),numConsFun);
numClust = zeros(1,numConsFun);

for c = 1:numConsFun
    params.HCL_clustMethod = consFun{c};
        
    [labelsCons(:,c), moreInfo] = pplk_runClusterer('HCL',WA, K, 1, params);
    numClust(c) = moreInfo{1}.numClust;
end

