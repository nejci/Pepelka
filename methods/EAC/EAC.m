function [labelsCons, numClust] = EAC(labelsEns,K,consFun)
% EAC
% [labelsCons, numClust] = EAC(labelsEns,K,consFun)
% Ensemble clustering by Evidence accumulation clustering (Fred & Jain, 2005)
%--------------------------------------------------------------------------
% INPUTS
%   labelsEns   (matrix)	clustering ensemble; each COLUMN corresponds
%							to one clustering (standard in Pepelka).
%
%   K			(scalar)	number of clusters in output consensus
%	                        clustering; leave empty [] to automatically
%	                        determine this number by lifetime criterion.
%
%   consFun     (string)    which clustering metod to use to cluster
%                           co-occurence matrix:
%                             'single' | 'SL'   - nearest distance (default)
%                             'complete' | 'CL' - furthest distance
%                             'average' | 'AL'  - unweighted average distance (UPGMA)
%                                             (also known as group average)
%                             'weighted'  - weighted average distance (WPGMA)
%                             'centroid'  - unweighted center of mass distance (UPGMC)
%                             'median'    - weighted center of mass distance (WPGMC)
%                             'ward'      - inner squared distance (min variance algorithm)
%               (cell)        combination of methods in cell string
%
% OUTPUTS
%   labelsCons    labels of consensus clustering(s)
%   numClust      number of clusters in consensus clustering
%--------------------------------------------------------------------------
% EXAMPLE
%     % Reproducing results from Fred & Jain, 2005; Table 3.
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
%     %% consensus using EAC
%     Kcons = []; % Change if you want to control the number of clusters in
%                 % the consensus partition.
% 
%     [labelsCons, numClust]= EAC(labelsEns,Kcons,'single');
% 
%     pplk_scatterPlot(data,labelsCons);
%     pplk_validExt(target,labelsCons,{'CA','NMI'})
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2014  Nejc Ilc
% Part of Pepelka package.
%------- VERSION ----------------------------------------------------------
% Version: 1.0
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
CO = computeCO(labelsEns,[],'labels','vec');

% Run linkage algorithm on CO vector/matrix to obtain consensus parititon.
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
    
    [labelsCons(:,c), moreInfo] = pplk_runClusterer('HCL',CO, K, 1, params);
    numClust(c) = moreInfo{1}.numClust;
end

% alternatively, one can run unwrapped function
% [labelsCons, moreInfo] = pplk_runClustererHCL(CO, K, consFun);