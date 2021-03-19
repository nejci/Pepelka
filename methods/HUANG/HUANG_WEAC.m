function [labelsCons numClust I_ncai] = HUANG_WEAC(labelsEns, K, consFun, betaOrINCAI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [labelsCons numClust I_ncai] = HUANG_WEAC(labelsEns, K, consFun, betaOrICNAI)
%
% INPUT
%     labelsEns - ensemble [N X ensembleSize]
%     K - desired number of clusters; leave [] for automatic determination
%     consFun - 'single', 'complete', 'average', 'ward' or cell string of
%     any combination of them
%     betaOrINCAI - parameter beta or vector of precomputed I_cnai values
%
% OUTPUT
%     labelsCons - consensus clustering
%     numClust - number of clusters in the labelsCons
%     I_ncai - vector of I_ncai values - if needed in another call
%--------------------------------------------------------------------------
% The source code of the WEAC method by Dong Huang.
% Version 1.0. May 14, 2014.
% Optimized for speed by Nejc Ilc, 10 June 2014
%
% If you use this code in your work, please cite the following paper:
%
% Dong Huang, Jian-Huang Lai, Chang-Dong Wang. 
% Combining Multiple Clusterings via Crowd Agreement Estimation and Multi-
% Granularity Link Analysis. Neurocomputing, in press, 2014.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('consFun','var') || isempty(consFun)
    consFun = 'single';
end
if ~exist('betaOrINCAI','var') || isempty(betaOrINCAI)
   betaOrINCAI = 2; 
end

if isscalar(betaOrINCAI)
    % Get the NCAI
    ncai = getNCAI(labelsEns);
    % The influence of the NCAI
    I_ncai = ncai.^betaOrINCAI;
else
    I_ncai = betaOrINCAI;
end


%  Get the weighted co-association matrix.
S = getWeightedMatrix(labelsEns, I_ncai);

% Run linkage algorithm on CO vector/matrix to obtain consensus parititon.
% If K is empty, number ob final clusters is automatically determined.
params = pplk_setParamsDefault();
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
    params.HCL_distance = [];
    
    [labelsCons(:,c), moreInfo] = pplk_runClusterer('HCL',S, K, 1, params);
    numClust(c) = moreInfo{1}.numClust;
end
