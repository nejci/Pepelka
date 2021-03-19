function value =indexDNg(data,labels,graph_type,options)
% INDEXDNG computes the Generalized Dunn's index on data using graph of
% graph_type. Best cluster partition maximizes the index value.
% value = INDEXDNG(data,labels,graph_type,options)
%--------------------------------------------------------------------------
% INPUTS:
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
%--------------------------------------------------------------------------
% OUTPUTS:
%   value        (scalar)	value of index
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2012  Nejc Ilc
%
%------- REFERENCE --------------------------------------------------------
% Pal, N. R., & Biswas, J. (1997). Cluster validation using graph theoretic
% concepts. Pattern Recognition, 30(6), 847-857.
%
%------- VERSION ----------------------------------------------------------
% Version: 1.1
% Last modified: 26-July-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

if (nargin < 2) 
	  error('indexDNg: Labels required!'); 
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

D = size(data,2);

% Construct graph on the data
[G,~,uniqueInd] = graph_create(data,labels,graph_type,options);

% if the removal of duplicates occured, update data and labels
if ~isempty(uniqueInd)
    data = data(uniqueInd,:);
    labels = labels(uniqueInd);
end

K=length(unique(labels));

%G=sqrt(G); %switch from squared euclidean to euclidean distances

% Compute diameter of each cluster.
% Diameter is the longest distance between any points x and y, where x and
% y belong to the same cluster.
diam=zeros(1,K);

% cluster centres - means 
center = zeros(K,D);

% diam(Ci) is maximum edge weight in the cluster Ci
% dist(Ci,Cj) is minimum distance between clusters centers (means?)

for k=1:K
	% select indices of nodes in cluster C_k
	C_ind=(labels==k);
	
	% pair-wise distances between all points in cluster C_k
	C_dist=G(C_ind,C_ind);	
	
	% find the highest value - longest distance in the cluster C_k
	diam(k)=max(max(C_dist));	
	
	% center is a mean of data points in cluster C_k
	center(k,:)=mean(data(C_ind,:),1);
end

max_diam=max(diam);

% distances between every pair of clusters i and j
dist_CiCj=Inf*ones(K,K);

for i=1:K-1
	for j=i+1:K
		% compute distance between C_i and C_j centres
		%dist_CiCj(i,j) = sqrt(dist_euclidean( center{i}, center{j}));
		dist_CiCj(i,j) = sqrt(sqdistance2( center(i,:), center(j,:)));
	end
end

value = min(min(dist_CiCj./max_diam));

end