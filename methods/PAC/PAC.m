function [labelsCons, numClust] = PAC(labelsEns,dim,K,consFun)
% [labelsCons, numClust] = PAC(labelsEns,dim,K,consFun)
% Probability accumulation clustering (Wang et al., 2009)
%--------------------------------------------------------------------------
% INPUTS
%   labelsEns   (matrix)	clustering ensemble; each COLUMN corresponds
%							to one clustering (standard in Pepelka).
%
%   dim         (scalar)    number of dimensions (features) in the data
%
%   K			(scalar)	number of clusters in output consensus
%	                        clustering; leave empty [] to automatically
%	                        determine this number by lifetime criterion.
%
%   consFun     (string)    which clustering metod to use to cluster
%                           co-occurence matrix:
%                             'single'    - nearest distance (default)
%                             'complete'  - furthest distance
%                             'average'   - unweighted average distance (UPGMA)
%                                             (also known as group average)
%                             'weighted'  - weighted average distance (WPGMA)
%                             'centroid'  - unweighted center of mass distance (UPGMC)
%                             'median'    - weighted center of mass distance (WPGMC)
%                             'ward'      - inner squared distance (min variance algorithm)
%               (cell)       combination of methods in cell string
%
% OUTPUTS
%   labelsCons    labels of consensus clustering(s)
%   numClust      number of clusters in consensus clustering
%--------------------------------------------------------------------------
% EXAMPLE
%     %% Load data
%     [data,target] = pplk_loadData('real\iris_orig');
% 
%     %% Generate ensemble
%     params = pplk_setParamsDefault();
%     params.KM_nRuns = 10;
%     params.dim = size(data,2);
%     ensSize = 50;
%     K = 3;
%     labelsEns = pplk_genEns(data,{'KM',ensSize,K,'fixed'},params);
% 
%     %% consensus using PAC
%     Kcons = []; % Change if you want to control the number of clusters in
%                 % the consensus partition.
% 
%     [labelsCons, numClust] = pplk_consEns(labelsEns,Kcons,'PAC-SL',params);
% 
%     pplk_scatterPlot(data,labelsCons);
%     pplk_validExt(target,labelsCons,{'CA','NMI'})
%
%
% Reference:
% X. Wang, C. Yang, and J. Zhou, �Clustering aggregation by probability
% accumulation,� Pattern Recognition, vol. 42, pp. 668�675, 2009.
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2014  Nejc Ilc
% Part of Pepelka package.
%------- VERSION ----------------------------------------------------------
% Version: 1.1
% Last modified: 23-June-2014 by Nejc Ilc
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

if ~exist('K','var')
    K=[]; % use lifetime cirteria to determine number of clusters
end

if ~exist('consFun','var') || isempty(consFun)
    consFun = 'single';
end

% Compute co-occurence matrix, return vector of pairs
PA = computePA(labelsEns,dim,'vec',[]);

% Run linkage algorithm on PA vector/matrix to obtain consensus parititon.
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
    
    [labelsCons(:,c), moreInfo] = pplk_runClusterer('HCL',PA, K, 1, params);
    numClust(c) = moreInfo{1}.numClust;
end

