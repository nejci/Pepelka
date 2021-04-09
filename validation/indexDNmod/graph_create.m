function [G,d,uniqueInd]=graph_create(data,labels,graph_type,options)
% GRAPH_CREATE  creates a graph of graph_type on data.
% Usage: [G,d,uniqueInd]=GRAPH_CREATE(data,labels,graph_type,options)
%--------------------------------------------------------------------------
% INPUTS:
%   data			(matrix)	matrix [n X d] with n d-dimensional samples
%	labels			(vector)	(for plotting only!) 
%                               array of non-negative integers determining
%								the labels of data samples
%                               use [] if there is no labels
%	graph_type		(scalar)	which type of graph to use. Select from:
%								'rng'		- relative neighborhood graph
%								'gabriel'	- Gabriel graph
%								'directedKnn' 
%								'mutualKnn'	
%								'symKnn'	- symmetric KNN graph
%								'EMST'		- Euclidean Minimum Spanning Tree graph
%								'epsilon'	- Epsilon graph
%	options.k - connectivity parameter of the kNN graph [default=3]
%	options.eps - threshold of the epsilon graph [default=auto]
%	options.show - if 0, plot of the graph is not shown [default=0]
%   options.graph_sqEucl - if 1, weights in graph are squared Euclidean
%                          distances, otherwise Euclidean (default).
%   options.graphSuper - super (global) graph on data'. It holds that data \subset data'
%   options.distMat - distance matrix (squared Euclidean distance)
%--------------------------------------------------------------------------
% OUTPUTS:
%   G - sparse matrix representing a graph
%   d - [n x n] symmetric matrix containing pair-wise distances between data points. 
%   uniqueInd - indices of unique data points if the duplicates removal
%               occurs; [] otherwise.
%--------------------------------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%--------------------------------------------------------------------------

if (nargin < 3) 
	  error(' GRAPH_CREATE: graph_type required!'); 
end

N = size(data,1);

if isempty(labels)
	labels=ones(N,1);
end

if ~exist('options','var')
	options = [];
end

if ~isfield(options,'k')
    options.k = min(size(data,1)-1,3);
end

if ~isfield(options,'eps')
    options.eps=[];
end

if ~isfield(options,'show')
    options.show=0;
end

if ~isfield(options,'graph_sqEucl')
    options.graph_sqEucl=0;
end

if ~isfield(options,'distMat')
    distMat = [];
else
	distMat = options.distMat;
end

if ~isfield(options,'graphSuper')
    options.graphSuper = [];
end



% Check data for duplicates and remove them.
[~,iA,~] = unique(data,'rows','first');
uniqueInd = [];

if length(iA) ~= N
    warning('graph_create:unique','There are duplicates in the data. Only unique data points will be considered.');
	if ~isempty(options.graphSuper)
		error('Super graph is defined but duplicates in sub data are found.');
	end
    % keep the same sequence of data - only delete duplicates
    iA = sort(iA);
    labels=labels(iA);
    data = data(iA,:);
    uniqueInd = iA;
end

if isempty(distMat)
	distMat = sqdistance2(data);
else
	if ~isempty(uniqueInd)
		distMat = distMat(iA,iA);
	end
end

% Construct graph on the data
% TODO: use distMat for other graph types as well (not only Gabriel)
switch(graph_type)
	case 'rng'
        [G,d] = graph_rng(data, 1e-10);
	
	case 'gabriel'
        %[G,d] = graph_gabriel(data, 1e-10, distMat, options.graphSuper);
        [G,d] = graph_gabriel(data, 1e-10, distMat);
	
	case 'EMST'
		[G,d] = graph_EMST(data,[]);
		
	% WARNING! Knn (directed, mutual, symm) and epsilon graph will probably 
    % consists of more than one component. Thus, there could be no 
    % connection between clusters. 
	case 'directedKnn'
		d = sqdistance2(data,data);
		G = graph_directedKnn(d,options.k,'dist');
	
	case 'epsilon'
        if options.graph_sqEucl == 0
            d = sqrt(sqdistance2(data,data));
        else
            d = sqdistance2(data,data);
        end		
		% if options.eps is [], estimation for "good" epsilon is calculated
		G = graph_epsilon(d,options.eps,'dist');
        
        
	case 'mutualKnn'
		d = sqdistance2(data,data);
		G = graph_mutualKnn(d,options.k,'dist');
	
	case 'symKnn'
		d = sqdistance2(data,data);
		G = graph_symmetricKnn(d,options.k,'dist');
			
	otherwise
		error(' GRAPH_CREATE: Wrong graph type!');
end

if options.graph_sqEucl == 0 && ~strcmpi(graph_type,'epsilon')
    G = sqrt(G); %switch from squared euclidean to euclidean distances
end

if (options.show==1)
	
	for L=1:size(labels,2)
		lab_uniq=unique(labels(:,L));
		K=length(lab_uniq);

		fig=figure();
		gplot(G,data,'-k');
		hold on;
		p_options.fig=fig;
		p_options.title=[graph_type, ' graph, labels=',num2str(L)];
		pplk_scatterPlot(data,labels(:,L),K,p_options);
		axis('equal');
		hold off;
	end
elseif (options.show==2)
	if strcmp(graph_type,'directedKnn')
		view(biograph(G,[],'ShowArrows','on','ShowWeights','on'));
	else
		view(biograph(G,[],'ShowArrows','off','ShowWeights','on'));
	end
end

end