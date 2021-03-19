function [labelsCons, K, time] = comusa(labelsEns, relaxation)
% [labelsCons, K, time] = COMUSA(labelsEns, relaxation)
%--------------------------------------------------------------------------
% COMUSA - COmbining Multiple clUsterings using Similarity grAph
% Consensus function for combining ensemble of clusterings.
%
% Based on Java implementation: http://www.cs.umb.edu/~smimarog/comusa/.
%--------------------------------------------------------------------------
% INPUTS
%   labelsEns		(matrix)	clustering ensemble; each COLUMN corresponds
%								to one clustering (standard in Pepelka).
%					(string)	name of the text file containing cluster
%								labels. One clustering per ROW or in a
%								binary format (each row corresponds to one
%								cluster)
%
%	relaxation		(scalar)	relaxation of maximum edge 
%								weight constraint. 0.5 means 50% relaxation 
%								and thus bigger final clusters.
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
% Mimaroglu, S. & Erdil, E. Combining multiple clusterings using similarity
% graph. Pattern Recognition, 2011, 44, 694 - 703.
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 24-April-2012 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <myName.mySurname@fri.uni-lj.si>
%==========================================================================



% determine whether labelsEns is filename or is it already a matrix.
if ischar(labelsEns)
	% PARSE
	% Each column corresponds to one data sample.
	% Objects are separated with comma.
	%--------------------------------------------------------------------------
	% labelsEns is a representation of multiple clusterings in a matrix form.
	labelsEns = csvread(labelsEns);
else
	labelsEns = labelsEns';
end

% In which format are labels organised (each row is clustering or is it
% cluster)?
% Convert to binary.
% Each row represents a cluster and each column shows cluster membership of
% an object in the cluster. Order of the clusters is not important.
if any(labelsEns > 1)
	[N,D]=size(labelsEns);
	PI = [];
	for ind=1:N
		m=length(unique(labelsEns(ind,:)>0));
		tmp=zeros(m,D);
		for d=1:D
			idx=labelsEns(ind,d);
			if idx > 0
				tmp(idx,d)=1;
			end
		end
		PI = [PI; tmp];
	end
else 
	PI = labelsEns;
end

PI = logical(PI);

% M = number of clusters
% N = number of data samples
[M,N] = size(PI);

tic
%--------------------------------------------------------------------------
% CONSTRUCT GRAPH = SM matrix
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SM is N x N co-association matrix. Cell (i,j) corresponds to the number
% of co-occurence of data sample i and j in the same cluster.
%--------------------------------------------------------------------------
SM = zeros(N); % TODO - could be sparse!

for i=1:N
	for j=i:N
		SM(i,j) = sum(PI(:,i) & PI(:,j));
	end
end
% copy upper triangular part of matrix to the lower one to make symmetric
% matrix 
SM = SM + triu(SM,1)';

% check for possible errors
if any(diag(SM) ~= sum(PI(:,1)))
	error('Wrong matrix SM!');
end

%--------------------------------------------------------------------------
% CALCULATE DEATTACHMENT
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Deattachment is calculated as df(i)/sw(i) and equals 1/attachment.
% df(i) = degree of freedom of a node (how many edges does node i have)
% sw(i) = sum of weights of those edges that are connected to node i
%--------------------------------------------------------------------------

df=sum(SM > 0,2)-1;
sw=sum(SM,2)-SM(1,1);
deattach = df./sw;

% remove self-loops
SM(linspace(1,numel(SM),length(SM))) = 0;

% if division by zero occurs
deattach(isnan(deattach)) = Inf;

%--------------------------------------------------------------------------
% MAIN PART - assigning to clusters
%--------------------------------------------------------------------------

% bit vector with ones where there is a visited node
visited = false(1,N);

% number of visited nodes
numVisited = 0;

% consensus clustering; labels vector and cell of clusters
labelsCons = zeros(N,1);
clusterID = 1;

% queue for pivots
Q = [];

% loop while there is a non-visited node
while numVisited ~= N
	
	% find a node with the lowest deattachment
	[~,pivot] = min(deattach);
	
	% add pivot to the queue and mark it as visited
	Q = [Q pivot];
	
	% mark node as visited
	visited(pivot) = 1; 
	numVisited = numVisited +1;
	
	while ~isempty(Q)
		
		% remove first element from Q
		v = Q(1); Q(1)=[];
		
		% to avoid selecting it again with min(deattach)
		deattach(v) = Inf;

		% add node to current cluster
		labelsCons(v) = clusterID;

		% expand cluster on neighbors
		%  -> find all non-visited neighbors of pivot node
		neighs = find((SM(:,v)>0) & (~visited)');
		
		%  -> loop over neighs of active pivot
		for nInd = 1:length(neighs)
			% one of the neighbors
			w = neighs(nInd);

			% check if w has other neighbor to which it is connected with
			% stronger weight.
			weight = SM(w,v);
			neighs2 = SM(:,w); % neighbor edges of w
			maxW = max(neighs2(~visited));

			% constraint
			if maxW > (weight + weight * relaxation)
				continue;
			end
			
			% add w to the queue
			Q = [Q w];
			
			% mark w as visited
			visited(w) = 1;
			numVisited = numVisited +1;			
		end
	end
	clusterID = clusterID + 1;
end
time = toc;
K = clusterID -1;



