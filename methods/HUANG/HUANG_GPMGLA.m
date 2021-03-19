function [labelsCons, numClust, I_ncai] = HUANG_GPMGLA(labelsEns, K, betaOrINCAI, alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [labelsCons numClust, I_ncai] = HUANG_GPMGLA(labelsEns, K, betaOrICNAI, alpha)
%
% INPUT
%     labelsEns - ensemble [N X ensembleSize]
%     K - desired number of clusters; leave [] for automatic determination
%     betaOrICNAI - parameter beta or vector of precomputed I_cnai values
%     alpha - parameter, default 0.5
%
% OUTPUT
%     labelsCons - consensus clustering
%     numClust - number of clusters in the labelsCons
%     I_ncai - vector of I_ncai values - if needed in another call
%--------------------------------------------------------------------------
% The source code of the GP-MGLA method by Dong Huang.
% Version 1.0. May 14, 2014.
% Optimized for speed by Nejc Ilc, 10 June 2014
%
% If you use this code in your work, please cite the following paper:
%
% Dong Huang, Jian-Huang Lai, Chang-Dong Wang.
% Combining Multiple Clusterings via Crowd Agreement Estimation and Multi-
% Granularity Link Analysis. Neurocomputing, in press, 2014.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimized for speed by Nejc Ilc, 5 June 2014

if ~exist('betaOrINCAI','var') || isempty(betaOrINCAI)
    betaOrINCAI = 2;
end
if ~exist('alpha','var') || isempty(alpha)
    alpha = 0.5;
end

if isscalar(betaOrINCAI)
    % Get the NCAI
    ncai = getNCAI(labelsEns);
    % The influence of the NCAI
    I_ncai = ncai.^betaOrINCAI;
else
    I_ncai = betaOrINCAI;
end

% Build bipartite graph
[bcs, baseClsSegs] = getAllSegs(labelsEns);
numClusters = size(baseClsSegs,1);
B = getBipartiteGraph(bcs, baseClsSegs, I_ncai, alpha);

if K > numClusters
    error('The cluster number in consensus clustering cannot exceeds the total number of clusters in the ensemble!!!');
else
    labelsCons = Tcut_for_clustering_ensemble(B,K);
end

numClust = length(unique(labelsCons));
