function [assignment_linkage,param_string]=cluster_single_linkage(S,k)
  
% Clusters the input similarity based on the single linkage as
% implemented in the statistics toolbox in matlab. 
%  
% INPUT:
% S = similarity matrix (called the Affinity Matrix in the paper)
% k = number of clusters
%
% OUTPUT
% assignment_linkage( 1,n ) vector of integers 1:k assigning the points to clusters
  
  if nargout ==2
	[assignment_linkage,param_string]=cluster_linkage(S,k,'single'); 
  else
	assignment_linkage=cluster_linkage(S,k,'single'); 
  end
  
  
  
