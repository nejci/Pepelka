function [CC,CCavg] = graph_localCC(G,nanMode)
% Computes local clustering coefficient
% G: adjancency matrix of a graph (directed or undirected)
% nanMode:  if 1, the CC of node with degree <= 1 is NaN and is not
%                 included in the average (CCavg).
%           if 0, the CC of node with degree <= 1 is 0 and IS
%                 included in the average (CCavg).
% CC: the local clustering coefficient of each node, 
%      where CC(i) = (num triangles connected to i) / (num triples centered on i)
% avgCC: average of all nodes except for those with degree <= 1.
%
% Ref: M. E. J. Newman, "The structure and function of complex networks"
% Written by Nejc Ilc, 2013
% Based on the code from Matlab Tools for Network Analysis
% Copyright (c) 2011, Massachusetts Institute of Technology.

if nargin < 2
   nanMode = []; 
end
if ~exist('nanMode','var') || isempty(nanMode)
   nanMode = 0; 
end

if nanMode
    degree_1_0 = NaN;
else
    degree_1_0 = 0;
end

n = length(G);
isdirected = ~isequal(G, G');

G = G > 0;  % convert to binary graph, no multiple edges
[deg,~,~] = graph_degrees(G,isdirected);
CC = zeros(n,1); % initialize clustering coefficient

% multiplication change in the clust coeff formula
coeff = 2;
% if directed graph
if isdirected
    coeff=1;
end

for i=1:n
    if deg(i) <= 1
        CC(i) = degree_1_0;
        continue;
    end    
    neigh = find(G(i,:)>0);    
    sub = G(neigh,neigh);  
    
    if isdirected
        edges_s = sum(sub(:));
    else
        sl = sum(diag(sub));
        edges_s = (sum(sub(:))-sl)/2+sl;
    end
    
    CC(i)=coeff*edges_s/deg(i)/(deg(i)-1);    
end

% Compute average CC
if nargout > 1
    CCv = CC(~isnan(CC));
    CCavg = mean(CCv);
end