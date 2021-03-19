function [labels, moreInfo]=pplk_clustererHCL(S,K,params)

% S:      - similarity matrix of size [N X N], where N is number of data points.
%         - vector of similarities of length N*(N-1)/2
%         - data matrix of size [N X noDimensions]; parameter distMethod
%           must be specified!
%
% K:      desired number of output clusters. If empty, K is automatically
%         determined using cluster lifetime criterion (Fred & Jain,2005)
%
% params.HCL_clustMethod: method for hierarchical clustering; choices are:
%       'single'    --- nearest distance (default)
%       'complete'  --- furthest distance
%       'average'   --- unweighted average distance (UPGMA) (also known as
%                       group average)
%       'weighted'  --- weighted average distance (WPGMA)
%       'centroid'  --- unweighted center of mass distance (UPGMC)
%       'median'    --- weighted center of mass distance (WPGMC)
%       'ward'      --- inner squared distance (min variance algorithm)
%
% params.HCL_distance: when S is a data matrix, this parameter must be specified
%       'euclidean'   - Euclidean distance (default)
%       'seuclidean'  - Standardized Euclidean distance. Each coordinate
%                       difference between rows in X is scaled by dividing
%                       by the corresponding element of the standard
%                       deviation S=NANSTD(X). To specify another value for
%                       S, use D=PDIST(X,'seuclidean',S).
%       'cityblock'   - City Block distance
%       'minkowski'   - Minkowski distance. The default exponent is 2. To
%                       specify a different exponent, use
%                       D = PDIST(X,'minkowski',P), where the exponent P is
%                       a scalar positive value.
%       'chebychev'   - Chebychev distance (maximum coordinate difference)
%       'mahalanobis' - Mahalanobis distance, using the sample covariance
%                       of X as computed by NANCOV. To compute the distance
%                       with a different covariance, use
%                       D =  PDIST(X,'mahalanobis',C), where the matrix C
%                       is symmetric and positive definite.
%       'cosine'      - One minus the cosine of the included angle
%                       between observations (treated as vectors)
%       'correlation' - One minus the sample linear correlation between
%                       observations (treated as sequences of values).
%       'spearman'    - One minus the sample Spearman's rank correlation
%                       between observations (treated as sequences of values).
%       'hamming'     - Hamming distance, percentage of coordinates
%                       that differ
%       'jaccard'     - One minus the Jaccard coefficient, the
%                       percentage of nonzero coordinates that differ
%       function      - A distance function specified using @, for
%                       example @DISTFUN.
%
%

Kdefined=1; % K is user-defined
if ~exist('K','var') || isempty(K)
    Kdefined = 0;
end

 
if ~exist('params','var') || isempty(params) || ~isfield(params,'HCL_clustMethod')
    clustMethod = 'single';
else
    clustMethod = params.HCL_clustMethod;
end

distMethod = 'euclidean';
if exist('params','var') && isfield(params,'HCL_distance')
    distMethod = params.HCL_distance;
end

% resolve synonyms
if strcmpi(clustMethod,'SL')
    clustMethod = 'single';
end
if strcmpi(clustMethod,'AL')
    clustMethod = 'average';
end
if strcmpi(clustMethod,'CL')
    clustMethod = 'complete';
end
if strcmpi(clustMethod,'WL')
    clustMethod = 'ward';
end


[N,Dim] = size(S);

% determine if S is a similarity vector or matrix or data matrix
Smode = 'data';
if ~exist('distMethod','var') || isempty(distMethod)
    Smode = 'simMat';
    if N == 1 || Dim == 1
        Smode = 'simVec';
    end
end

if strcmp(Smode,'simMat') && N ~= Dim
    error('Matrix S must be either a square similarity matrix or data matrix - in the latter case, you must specify parameter distMethod!');
end

%==========================================================================
tic


if strcmp(Smode,'data')
    htree=linkage(S,clustMethod,distMethod);
    
else
    
    MAX = max(S(:));
    MIN = min(S(:));
    
    % Normalize similarities on [0,1].
    if MAX > 1 || MIN < 0
        S = (S-MIN) ./ (MAX-MIN);
    end
    
    % convert similarities to distances
    D = 1 - S;
    
    if strcmp(Smode,'simMat')
        % reformat distance matrix to a row vector of length m(m–1)/2,
        % corresponding to pairs of observations in a matrix X with m rows
        % Distances arranged in the order (2,1), (3,1), ..., (m,1), (3,2), ...,
        % (m,2), ..., (m,m–1))
        V = nchoose2(1:N);
        I = V(:,1)';
        J = V(:,2)';
        ind=(I-1)*N+J;
        
        D = D(ind);
    end
    htree=linkage(D,clustMethod);
end

% Determine output labels
% If desired number of clusters is not defined, automaticaly find it using
% cluster lifetime criterion
if Kdefined
    labels=cluster(htree,'maxclust',K);
else
    % lifetime criterion
    %dendrogram(htree,0);
    lifetime= htree(2:end,3)-htree(1:end-1,3);
    [~, lifeT_i] = max(lifetime);
    K = length(lifetime)+2 - lifeT_i;   
    
    %  break tree into clusters by cutoff
    c = (htree(lifeT_i+1,3) + htree(lifeT_i,3)) /2;
    labels = cluster(htree,'cutoff',c,'criterion','distance');

    if K ~= length(unique(labels))
       warning('pplk:HCL','Something wrong with a tree-cut.'); 
    end    
end
time=toc;

moreInfo.time=time;
moreInfo.numClust = K;