% EMST: Finds the Euclidean Minimum Spanning Tree graph among points in P dimensions.
%
%     Usage: [G,dist] = graph_EMST(data)
%
%           data =    [n x p] matrix of point coordinates.
%           noplot =  optional boolean flag indicating that plot should be 
%                       suppressed [default = 0].  Plot of the Gabriel graph
%                       is produced only for p=2.
%           tol =     optional tolerance for squared distance from a third point 
%                       to the minimum A-B circle (i.e., that circle on whose 
%                       circumference A & B are at opposite points) [default = 1e-6].
%           -------------------------------------------------------------------------
%			G =		  [n x n] sparse matrix that represents a graph. Nonzero
%					  entries in matrix G represent the weights of the edges.
%           connect = [n x n] boolean adjacency matrix.
%           dist =    corresponding edge lengths (Euclidean distances);
%                       non-connected edge distances are given as zero.
%
% Implementation inspired by
% http://en.wikipedia.org/wiki/Euclidean_minimum_spanning_tree

function [G,dist,connect] = graph_EMST(data,options)
  
if nargin < 2
	options.show=[];
end

if isempty(options)
	options.show=0;
end
[N,D]=size(data);

% compute Delaunay triangulation
% It is very slow for higher dimensional cases (D>3), therefore it is not
% supported (yet).

% 2-D and 3-D case
if D < 4
	DT = DelaunayTri(data);
else
	%DT = delaunayn(data);
    error('EMST graph for data with more than 3 dimensions is not supported.');
end

E = edges(DT);
nEdges=size(E,1);
dist = sqdistance2(data);

%create sparse matrix - graph representation
C=sparse(E(:,1),E(:,2),ones(nEdges,1),N,N,nEdges);
connect=C;
G= C+C';
% weight graph by distances
G=G.*dist;

% graph is undirected, so we can forget half of edges
G=tril(G);

% compute minimal spanning tree on graph
[G] = graphminspantree(G);

if (options.show==1)
	fig=figure();
	gplot(G,data,'-k');
	hold on;
	plot(data(:,1),data(:,2),'.k');
	axis('equal');
end

