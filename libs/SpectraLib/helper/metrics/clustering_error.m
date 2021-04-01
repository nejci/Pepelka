function [min_perm_error,permutation] = clustering_error( true_assign, assign )
% function [min_perm_error,permutation] = clustering_error( true_assign, assign )
%
% Computes the clustering error in assign w.r.t. true_assign. The assign
% MUST be of same size 
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
  confusion=confusion/n; 
  BIPARTITE_METHOD=1; 
  if BIPARTITE_METHOD
	[max_val, matching] =max_bipartite_matching(confusion); 
	
	min_perm_error=1-max_val;
	if nargout>1
	  permutation=matching;
	end
  else 

	
	% The error is obv the sum of non-diagonal entries. However to take
	% into account the fact the cluster '1' maybe be clustered  as cluster
	% '2' the error needs to minimized over all the permutations of
	% assign.

	
	
	%generate all the assignments of assigned_k to true_k 
	all_combinations=combnk(1:k_all,k_true); 
	n_c_k=nchoosek(k_all,k_true);
	k_fact=factorial(k_true);
	num=n_c_k*k_fact;  
	all_permutations=zeros(num,k_true);
	for i=1:n_c_k
      start_lim=(i-1)*k_fact+1;
      all_permutations(start_lim:start_lim+k_fact-1,:)=perms(all_combinations(i,:));
	end
	
	max_correct_classification=0;
	for i=1:num
	  correct_classification=0;
	  for j=1:k_true
		correct_classification = correct_classification+ confusion(j, ...
												  all_permutations(i,j));
	  end
	  if correct_classification>max_correct_classification
		max_correct_classification=correct_classification;
		perm=all_permutations(i,:);
	  end
	end
	DIFF=max_correct_classification-max_bipartite_matching(confusion);
	if DIFF>1e-8
	  disp('BAD DIFF'); 
	end
	
	min_perm_error=1-max_correct_classification;
	if nargout>1
	  permutation=perm;
	end
	
  end
  
