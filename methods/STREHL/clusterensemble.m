% function cl = clusterensemble(cls,k,method,weights)
%
% DESCRIPTION
%   - performs the supra-consensus function for CLUSTER ENSEMBLES
%     to combine multiple clustering label vectors e.g., cls = [cl1; cl2;...]
%   - returns a row vector of the combined clustering labels
%   - the following consensus functions are computed 
%     - Cluster-based Similiarty Partitioning Algorithm (CSPA)
%     - HyperGraph Partitioning Algorithm (HGPA)
%     - Meta-CLustering Algorithm (MCLA)
%     and the one with the maximum average normalized mutual information
%     is returned
% ARGUMENTS
%   cls   - matrix of labelings, one labeling per row (n x p)
%           entries must be integers from 1 to k1 (row 1), 1 to k2 (row 2),...
%           use NaN as an entry for unknown/missing labels
%   k     - 1,2,3,... maximum number of clusters in the combined clustering
%           (optional, default max(max(cls))
%	method - CSPA|HGPA|MCLA or all of them if method=[] or non-existent
%
% EXAMPLE
%   clusterensemble; 
%   clusterensemble([ones(3,20) 2*ones(3,30); ones(1,10) 2*ones(1,40)])
% REFERENCE
%   please refer to the following paper if you use CLUSTER ENSEMBLES
%     A. Strehl and J. Ghosh, "Cluster Ensembles - A Knowledge Reuse
%     Framework for Combining Partitionings", Proc. of 18th National
%     Conference on Artificial Intelligence (AAAI 2002), July 2002,
%     Edmonton, Canada
% RELEASE
%   version 1.0, 2002/04/20, tested on Matlab 5.2.0.3084, LNX86
%   available from http://www.strehl.com
%   license granted for research use ONLY (see README)
%   copyright (c) 1998-2002 by Alexander Strehl

function cl = clusterensemble(cls,k,method,weights)

if ~exist('cls','var')
    disp('clusterensemble-warning: no arguments - displaying illustrative example:');
    cls = [1 1 1 2 2 3 3;2 2 2 3 3 1 1;1 1 2 2 3 3 3;1 2 NaN 1 2 NaN NaN];
    disp('clusterensemble-advice: type "help clusterensemble" for information about usage');
    disp(' ');
end

if ~exist('method','var') || isempty(method)
	if size(cls,2)>1000,
	  workfcts = {'hgpa', 'mcla'};
	  disp('STREHL clusterensemble.m - WARNING: using only hgpa and mcla because of speed.');
	else
	  workfcts = {'cspa', 'hgpa', 'mcla'};
	end
else
	workfcts={lower(method)};
end

if ~exist('weights','var')
    weights = [];
end
if ~exist('k','var')
    k = [];
end

numFcts = length(workfcts);
cl = zeros(numFcts,size(cls,2));
for i = 1:numFcts
   workfct = workfcts{i};

   cl(i,:) = feval(workfct,cls,k,weights);
   
   if length(workfcts) > 1
       q(i) = ceevalmutual(cls,cl(i,:));
   end
   %disp(['clusterensemble: ' workfct ' at ' num2str(q(i))]);
end

%supra-consensus function, if more than 1 methods are chosen
if length(workfcts) > 1
	[qual, best] = max(q);
	disp(['Strehl - Supra consensus function chose: ', workfcts{best},', NMI: ',num2str(qual)]);
	cl = cl(best,:);
end
