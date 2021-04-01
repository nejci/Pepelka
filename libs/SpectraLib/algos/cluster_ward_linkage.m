function [assignment_ward,param_string]=cluster_ward_linkage(S,k)
  
% Clusters the input similarity based on the ward linkage as
% implemented in the statistics toolbox in matlab. 
%  
% INPUT:
% S = similarity matrix (called the Affinity Matrix in the paper)
% k = number of clusters
%
% OUTPUT
% assignment_ward( 1,n ) vector of integers 1:k assigning the points to clusters
  
  if nargout ==2
	[assignment_ward ,param_string]=cluster_linkage(S,k,'ward'); 
  else
	assignment_ward=cluster_linkage(S,k,'ward'); 
  end
  
  
  
