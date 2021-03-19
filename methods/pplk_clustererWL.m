function [labels, moreInfo]=pplk_clustererWL(data,K,params)

% dist can be:
%         'euclidean'   - Euclidean distance
%         'seuclidean'  - Standardized Euclidean distance, each coordinate
%                         in the sum of squares is inverse weighted by the
%                         sample variance of that coordinate
%         'cityblock'   - City Block distance
%         'mahalanobis' - Mahalanobis distance
%         'minkowski'   - Minkowski distance with exponent 2
%         'cosine'      - One minus the cosine of the included angle
%                         between observations (treated as vectors)
%         'correlation' - One minus the sample linear correlation between
%                         observations (treated as sequences of values).
%         'spearman'    - One minus the sample Spearman's rank correlation
%                         between observations (treated as sequences of values).
%         'hamming'     - Hamming distance, percentage of coordinates
%                         that differ
%         'jaccard'     - One minus the Jaccard coefficient, the
%                         percentage of nonzero coordinates that differ
%         'chebychev'   - Chebychev distance (maximum coordinate difference)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modification: 14.8.2013 
% (C) Pepelka Package, Nejc Ilc (nejc.ilc@fri.uni-lj.si)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('params','var') && isstruct(params) && isfield(params,'WL_distance')
    dist=params.WL_distance;
else
	dist='euclidean';
end

tic;
D=pdist(data,dist);
htree=linkage(D,'ward'); 
labels=cluster(htree,K); 
time=toc;

moreInfo.time=time;
end