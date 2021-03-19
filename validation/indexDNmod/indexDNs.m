function value = indexDNs(data,labels,graph_type,options,penalFun,penalFunP,penalFunLogMid)
% value = INDEXDNS(data,labels,graph_type,options)
%--------------------------------------------------------------------------
% Modified Dunn's cluster validity index using shortest paths.
% Best cluster partition maximizes the index value.
%--------------------------------------------------------------------------
% INPUTS
%   data		(matrix)	matrix [n X d] with n d-dimensional samples
%
%	labels		(vector)	array of non-negative integers determining
%								the labels of data samples
%
%	graph_type	(scalar)	which type of graph to use. Select from:
%                             'rng'		    - relative neighborhood graph
%                             'gabriel'	    - Gabriel graph
%                             'directedKnn' - directed kNN graph
%                             'mutualKnn'	- mutual kNN graph
%                             'symKnn'	    - symmetric KNN graph
%                             'EMST'		- Euclidean Minimum Spanning
%                                             Tree graph
%                             'epsilon'	    - Epsilon graph
%
%	options.k    : connectivity parameter of the kNN graph [default=3]
%	options.eps  : threshold of the epsilon graph [default=auto]
%	options.show : if 0, plot of the graph is not shown [default=0]
%   options.graph_sqEucl - if 1, weights in graph are squared Euclidean
%                          distances, otherwise Euclidean.
%   penalFun        penal function for increasing K (only for DNS)
%                   ['none', 'reciproc', 'exponent', default: 'logistic']
%   penalFunP       strength of penals [0,1], 0: no penal, 1: full penal
%   penalFunLogMid  mid-point for logistic function, default ceil(sqrt(N)/2)
%--------------------------------------------------------------------------
% OUTPUTS:
%   value        (scalar)	value of index
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2016,  Nejc Ilc
%
%------- REFERENCE --------------------------------------------------------
% Ilc, N. (2012). Modified Dunn�s cluster validity index based on graph
% theory. Przeglad Elektrotechniczny (Electrical Review), (2), 126-131.
%
%------- VERSION ----------------------------------------------------------
% Version: 2.0
% Last modified: 18-Feb-2016 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================


if (nargin < 2)
    error('indexDNs: Labels required!');
end

if (nargin < 3)
    graph_type = [];
end
if  ~exist('graph_type','var') || isempty(graph_type)
    graph_type='gabriel';
end

if ~exist('options','var') || isempty(options)
    options.show = 0;
    options.eps = [];
    options.k = 3;
    options.graph_sqEucl = 0;
end

% Construct graph on the data
[G,~,uniqueInd] = graph_create(data,labels,graph_type,options);
isDirected = strcmpi(graph_type,'directedKnn');

% if the removal of duplicates occured, update labels
if ~isempty(uniqueInd)
    data = data(uniqueInd,:);
    labels = labels(uniqueInd);
end
% fix inconsistent labels
[~,~,labels] = unique(labels);
K = max(labels);
N = length(labels);
assert(N==size(data,1),'Length of labels and size of data matrix incosistent.');


% Find all (shortest) distances between data points (nodes in the [un]directed graph G)
[dist_global] = graphallshortestpaths(G,'Directed',isDirected);

% Compute diameter of each cluster.
% Diameter is the longest distance between any points x and y, where x and
% y belong to the same cluster.
diam = zeros(1,K);

for k = 1:K
    %select indices of data points in cluster C_k
    C_ind = (labels==k);
    
    % check for singletons and two-points cluster and bypass graph stuff
    C_size = sum(C_ind);
    if C_size == 1
        diam(k) = 0;
    elseif C_size == 2
        dataSub = data(C_ind,:);
        diam(k) = sqrt(sum((dataSub(1,:)-dataSub(2,:)).^2));
    else
        % create new graph only with vertices from the current cluster
        Gsub = graph_create(data(C_ind,:),[],graph_type,options);
        
        % Find all (shortest) distances between data points (nodes in the [un]directed graph G)
        % Restrict search only on edges inside cluster = sub graph
        C_dist = graphallshortestpaths(Gsub,'Directed',isDirected);
        
        %pair-wise distances between all points in cluster C_k
        %some entries of C_dist matrix can contain Inf value - it means that
        %there is no path between two points in cluster. These Inf values can
        %be ignored - we can set them to NaN.
        C_dist(C_dist==Inf)=NaN;
        if sum(sum(isnan(C_dist)))>0
            warning('INDEX_DNS:unconnected','Some points in the graph are not connected to the rest of the cluster! Decrease tolerance when creating graph.');
        end
        diam(k) = max(max(C_dist));
    end
end

max_diam = max(diam);

% distances between every pair of clusters i and j
dist_CiCj=Inf*ones(K,K);
for i=1:K-1
    for j=i+1:K
        % distance between clusters C_i and C_j is minimum distance between
        % any two points x and y, where x \in C_i and y \in C_j
        Ci_members = (labels==i);
        Cj_members = (labels==j);
        
        % select distances between points in C_i and C_j clusters.
        CiCj = dist_global(Ci_members,Cj_members);
        
        % compute minimum distance between C_i and C_j
        dist_CiCj(i,j) = min(min(CiCj));
        
        % if C_i and C_j are not connected,
        if isinf(dist_CiCj(i,j))
            warning('INDEX_DNS:unconnected','Clusters are not connected!');
        end
    end
end
min_dist_CiCj = min(dist_CiCj(:));

% penalities on K
% PenalFun arguments check
if ~exist('penalFun','var') || isempty(penalFun)
    penalFun = 'logistic';
end
if ~exist('penalFunP','var') || isempty(penalFunP)
    penalFunP = 1;
end
if ~exist('penalFunLogMid','var') || isempty(penalFunLogMid)
    penalFunLogMid = ceil(sqrt(N)/2);
end

p = penalFunP;
switch penalFun
    case 'none'
        penal = 1;
    case 'reciproc'
        penal = 1/(K^p);
    case 'exponent'
        penal = exp(-K*p);
    case 'logistic'
        K0 = penalFunLogMid; % sigmoid midpoint
        penal = 1 / (1 + p*exp(p*(K-K0)));
    otherwise
        error('Wrong penalize function.');
end

% Value of the DNS index
value = min_dist_CiCj / max_diam * penal;

end