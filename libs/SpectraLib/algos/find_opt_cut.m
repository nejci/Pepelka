function [min_cut_val, orig_indices_first_cluster]=find_opt_cut(origA, indices, origS,scaling_option,optimization_option) 
% function [min_cut_val, orig_indices]=find_cut_opt(origA, indices, origS) 
% helper function for Cluster_kvv  
  n=length(indices); 
  B=origA(indices,indices); 
  Browsum=sum(B,2);
  switch (scaling_option) 
   case {'add_diagonal'}
	B=B+diag(1-sum(B,2)); % augment the diagonal entries so as to % make the row sum 1. 
   case {'scale_row'}
	for i=1:n
	  B(i,:)=B(i,:)/Browsum(i); 
	end
   otherwise
	warning('******** WRONG scaling_option in cluter_kvv ********'); 
  end
  
  [vv lambda]=myeigs(B,2); 
  v2=vv(:,2); 
  for i=1:n
	dot_prod(i)=dot(v2,B(:,i)');
  end
  
  [dummy, local_ordering]=sort(dot_prod); 
  orig_ordered_indices=indices(local_ordering); 
  [min_cut_point,min_cut_val]= find_opt_cut_ordered(origS(orig_ordered_indices,orig_ordered_indices),optimization_option);
  
  first_cluster_local_indices=local_ordering(1:min_cut_point); 
  orig_indices_first_cluster=indices(first_cluster_local_indices); 
  
  
