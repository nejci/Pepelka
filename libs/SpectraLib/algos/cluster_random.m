function [cluster_assignment,param_string]=cluster_random(S,k)

% function cluster_assignment=cluster_random(S,k,nmin_mstcluster)
%  
% Random assignment 
%
% INPUT
%   S = similarity matrix
%   k = number of clusters
%
% OUTPUT
%   cluster_assignment  (1, Ndata) assignment of points to clusters 0:k
  
  n=length(S); 
  ca=rand(n,1);
  cluster_assignment=floor(1+k*ca);
  if nargout >1 
	param_string='none';
  end
  