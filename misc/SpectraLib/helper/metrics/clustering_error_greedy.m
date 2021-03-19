function [min_perm_error,permutation] = clustering_error_greedy( true_assign, assign )

%function distance = compare_clusterings( true_assign, assign )
%
% Computes the clustering error in assign w.r.t. true_assign. The assign
% MUST be of same size. (Uses the greedy approach and not minimum) 
%  
%
%
%true_assign, assign(1,n) = vectors of integers from 1:k that represent the
%                      assignment to clusters. must have the same length
%
% Note that this would just ignore all the points which have teh
% true_assign as zero. These are treated as outliers.


  %compute the confusion matrix and some of the parameters that are
  %used. 
  [n, k_true, k_all,confusion]=compute_confusion(true_assign,assign);

  
  min_perm_error=1;
  % Greedy approach 
  for i=1:k_true
	[val_array indices] = max(confusion); 
	[maxval maxcol]=max(val_array); 
	ki=indices(maxcol); 
	kalli=maxcol; 
	permutation(ki)=kalli;
	confusion(ki,:)=-1; 
	confusion(:,kalli)=-1; 
	min_perm_error=min_perm_error-maxval; 
  end
  
  % The error is obv the sum of non-diagonal entries. However to take
  % into account the fact the cluster '1' maybe be clustered  as cluster
  % '2' the error needs to minimized over all the permutations of
  % assign.

