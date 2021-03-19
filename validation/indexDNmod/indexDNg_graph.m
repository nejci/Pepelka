function value =indexDNg_graph(G,labels,data)
% INDEXDNG_GRAPH computes the Generalized Dunn's index on already built 
% graph G. Best cluster partition maximizes the index value.
% value = INDEXDNG_GRAPH(G,labels,data)
%--------------------------------------------------------------------------
% INPUTS:
%   G       (matrix)	graph on data points; output of function 
%                       graph_create()
%
%	labels	(vector)	array of non-negative integers determining
%						the labels of data samples
%
%   data	(matrix)	matrix [n X d] with n d-dimensional  samples,
%                       on which the graph G was constructed
%--------------------------------------------------------------------------
% OUTPUTS:
%   value	(scalar)	value of index
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013  Nejc Ilc
%
%------- REFERENCE --------------------------------------------------------
% Pal, N. R., & Biswas, J. (1997). Cluster validation using graph theoretic
% concepts. Pattern Recognition, 30(6), 847-857.
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 11-June-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================

if (nargin < 2) 
	  error('indexDNg: Labels required!'); 
end

K=max(labels);

% Compute diameter of each cluster.
% Diameter is the longest distance between any points x and y, where x and
% y belong to the same cluster.
diam=zeros(1,K);

% cluster centres - means 
center=cell(1,K);

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
	center{k}=mean(data(C_ind,:),1);
end

max_diam=max(diam);

% distances between every pair of clusters i and j
dist_CiCj=Inf*ones(K,K);

for i=1:K-1
	for j=i+1:K
		% compute distance between C_i and C_j centres
		dist_CiCj(i,j) = sqrt(sqdistance2( center{i}, center{j}));
		%dist_CiCj(i,j) = sqdistance2( center{i}, center{j});
	end
end

value = min(min(dist_CiCj./max_diam));

end