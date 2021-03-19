function [labelsCons, K, time] = diclens(labelsEns)

% [labelsCons, K, time] = DICLENS(labelsEns)
%--------------------------------------------------------------------------
% DICLENS: Divisive Clustering Ensemble with Automatic Cluster Number
% Consensus function for combining ensemble of clusterings.
%
% Based on Java implementation: http://www.cs.umb.edu/~smimarog/diclens/.
%--------------------------------------------------------------------------
% INPUTS
%   labelsEns		(matrix)	clustering ensemble; each COLUMN corresponds
%								to one clustering (standard in Pepelka).
%					(string)	name of the text file containing cluster
%								labels. One clustering per ROW or in a
%								binary format (each row corresponds to one
%								cluster)
%--------------------------------------------------------------------------
% OUTPUTS:
%   labelsCons		(vector)	ensemble clustering result - labels of
%								consensus clustering in a column vector
%	K				(scalar)	number of clusters in consensus clustering
%	time			(scalar)	execution time in seconds
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2012  Nejc Ilc
% Part of Pepelka package.
%
%------- REFERENCE --------------------------------------------------------
% Mimaroglu, S. & Aksehirli, E. 
% DICLENS: Divisive Clustering Ensemble with Automatic Cluster Number. 
% Computational Biology and Bioinformatics, IEEE/ACM Transactions on, 
% 2012, 9, 408-420
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 2-July-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <myName.mySurname@fri.uni-lj.si>
%==========================================================================

if ischar(labelsEns)
%--------------------------------------------------------------------------
% PARSE
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Each column corresponds to one data sample.
% Each row represents a clustering and each column shows cluster membership of
% an object in the cluster. Order of the rows is not important.
% Objects are separated with commas.
% Unknown values are given as 'N'.
%--------------------------------------------------------------------------
	file = fileread(labelsEns);
	file = strrep(file,'N','NaN');

	% PI is an representation of multiple clusterings in a matrix form.
	PI = str2num(file);
else
	PI = labelsEns';
end

tic;

% M = number of clusterings
% N = number of data samples
[M,N] = size(PI);

nClustL = max(PI,[],2);
nClust = sum(nClustL);

lookUpTable = zeros(nClust,3);
lookUpTable(:,1)=1:nClust;
ind=[0; cumsum(nClustL)];

for cInd=1:M
	lookUpTable(ind(cInd)+1:ind(cInd+1),2) = cInd;
	lookUpTable(ind(cInd)+1:ind(cInd+1),3) = 1:nClustL(cInd);
end

P = []; % clusters
% Construct binary representation of PI and save it to P.
% Each cluster is a row with ones on indeces of containing data points.
for p=1:M
	Ptmp = false(nClustL(p),N);
	
	for d=1:N
		ind=PI(p,d);
		if ~isnan(ind)
			Ptmp(ind,d)=1;
		end
	end
	P = [P; Ptmp];
end

%--------------------------------------------------------------------------
% Construct graph
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Graph is represented by sparse matrix G, where each entry determine an
% edge and weight at the same time. Weights of the edges are all nonzero 
% entries in the lower triangle of the N-by-N sparse matrix G.
% G(1,2) = 0.5 means that there is edge
% between node 1 and node 2 and is weighted by 0.5.
%--------------------------------------------------------------------------
G = sparse(nClust,nClust);
for ci=1:nClust
	for cj=ci+1:nClust
		% we switch indeces to produce lower triangular matrix G
		G(cj,ci) = computeECS(P(ci,:),P(cj,:),P,M);
	end
end

% cost of edges and weights have inverse relationship:
% Ci = max(W) - Wi + 1
maxW = max(max(G)) + 1;
%Gtrans = maxW - G + 1;

Gtrans = sparse(nClust,nClust);
for ci=1:nClust
	for cj=ci+1:nClust
		% transform weights to cost
		Gtrans(cj,ci) = maxW - G(cj,ci);
	end
end

[SMST, pred] = graphminspantree(Gtrans, 'Method', 'Prim');
%view(biograph(SMST,cellstr(int2str(lookUpTable(:,2:3))),'ShowArrows','off','ShowWeights','on'))

% sort edges in ascending order, from minimal to maximal cost
[Eind_i, Eind_j,E_val] = find(SMST);
E_ind =[Eind_i, Eind_j];
[E_val,pInd] = sort(E_val,'ascend');
E_ind = E_ind(pInd,:);

% for debug
%[lookUpTable(E_ind(:,1),[2:3]) lookUpTable(E_ind(:,2),[2:3]) maxW - E_val];

E_len=length(E_ind);

metaClusterGraph = sparse(nClust, nClust);
avgICS = zeros(1,E_len);
avgECS = zeros(1,E_len);
compHist = zeros(E_len,nClust); % history of component's labels

for edgeInd = 1 : E_len

	% poiscemo povezane komponente v metaClusterGraph
	[numComp, comp] = graphconncomp(metaClusterGraph, 'Directed', false);
	compHist(numComp,:) = comp;
		
	finalClusters = majorityVoting(comp, numComp, P);
	
	[nFin,N] = size(finalClusters);
	
	tmp=0;
	for cInd = 1:nFin
		tmp = tmp + computeICS(finalClusters(cInd,:), P);
	end
	avgICS(numComp) = tmp / nFin;
	
	tmp=0;
	for cInd1 = 1:nFin
		for cInd2 = cInd1+1:nFin
			tmp = tmp + computeECS(finalClusters(cInd1,:), finalClusters(cInd2,:), P, M);
		end
	end
	avgECS(numComp) = tmp / (nFin * (nFin-1) /2);
	
	metaClusterGraph(E_ind(edgeInd,1),E_ind(edgeInd,2)) = E_val(edgeInd); 
end	

% normalize ICS and ECS on interval [0,1]
avgICS = normalize_minmax(avgICS);
avgECS = normalize_minmax(avgECS);

% quality function is difference between normalized ICS and ECS
differ = avgICS - avgECS;

% the solution is the partition of metacluster graph with the highest differ
% value
[v,maxInd] = max(differ);

finalClusters = majorityVoting(compHist(maxInd,:), maxInd, P);
[v,labelsCons]=max(finalClusters,[],1);

time = toc;
labelsCons = labelsCons';
K = size(finalClusters,1);



