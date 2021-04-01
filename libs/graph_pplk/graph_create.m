function [G,d,data,labels]=graph_create(data,labels,graph_type,options)

% data = [n x d] matrix of n data points in d-dimensional space
% labels = (for plotting only!) optional numeric labels of nodes (1:K, where K is number of
% groups), use [] if there is no labels
% graph_type =	'rng' | 'gabriel' | 'directedKnn' | 'mutualKnn' | 'symKnn' |
%				'EMST' | 'epsilon'
% options.k - connectivity parameter of the kNN graph
% options.eps - threshold of the epsilon graph
% options.show - if 0, plot of the graph is not shown
%
% G - sparse matrix representing a graph
% d - [n x n] symmetric matrix containing pairwise distances between data points. 
% data - may be different from input data due to the duplicates removal
% labels - may be different from input labels due to the duplicates removal

if (nargin < 3) 
	  error(' GRAPH_CREATE: graph_type required!'); 
end
if (nargin < 4) 
	  options = [];
end
if isempty(options)
	options.eps = [];
	options.k=3;
	options.show=0;
end


% remove duplicates
%[data,m,n]=unique(data,'rows');
%labels=labels(m,:);

[N,dim] = size(data);
if isempty(labels)
	labels = ones(N,1);
end


% Construct graph on the data
switch(graph_type)
	case 'rng'
		[G,d] = graph_rng(data,1e-10);
	
	case 'gabriel'
		[G,d] = graph_gabriel(data,1e-10);
        
	case 'EMST'
        if dim > 3
            error('EMST graph for data with more than 3 dimensions is not supported.');
        end
		[G,d] = graph_EMST(data,[]);
		
	% WARNING! Knn (directed, mutual, symm) and epsilon graph will probably consists of more than one component,
	% thus there could be no connection between clusters. 
	case 'directedKnn'
		d = sqdistance2(data);
		G = graph_directedKnn(d,options.k,'dist');
	
	case 'epsilon'
		d = sqdistance2(data);
		% if options.eps is [], estimation for "good" epsilon is calculated
		G = graph_epsilon(d,options.eps,'dist');
	
	case 'mutualKnn'
		d = sqdistance2(data);
		G = graph_mutualKnn(d,options.k,'dist');
	
	case 'symKnn'
		d = sqdistance2(data);
		G = graph_symmetricKnn(d,options.k,'dist');
			
	otherwise
		error(' GRAPH_CREATE: Wrong graph type!');
end

%G=sqrt(G); %switch from squared euclidean to euclidean distances

if (options.show==1)
	
	for L=1:size(labels,2)
		lab_uniq=unique(labels(:,L));
		K=length(lab_uniq);

		%p_options.fig=fig;
		p_options.title=[graph_type, 'labels=',num2str(L)];
		fig = pplk_scatterPlot(data,labels(:,L),K,p_options);
        figure(fig);
        hold on;
        gplot(G,data,'-k');
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