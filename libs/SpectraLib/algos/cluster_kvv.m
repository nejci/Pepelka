function [assignment_kvv ,param_string]=cluster_kvv(S,k,scaling_option, optimization_option) 

% function [assignment_kvv param_string]=cluster_kvv(S,k) 
% 
% Cluster KVV : Based on the algorithms mentioned in Kannan Vempala and
% Vetta in "On Clustering : Good, Bad and Spectral"
  
  P=S_to_P(S); 
  n=length(S); 
  
  min_cut_val_arr = -ones(1,n); 
  cluster_assignment=ones(1,n); % put all the points in one cluster
                                % first. 
 
  nc=1; % number of clusters till now
  
  MIN_CLUST_SIZE=4;
  while (nc  < k )
	% Have nc clusters already assigned and want to split one into two to
    % have nc+1 clusters. 
	
	%ensure that all min_cuts are computed. 
	
	for i=1:nc
	  if min_cut_val_arr(i)==-1 
		indices=find(cluster_assignment==i); 
		if length(indices) > MIN_CLUST_SIZE
		  [min_cut_val, orig_indices_first_cluster]=find_opt_cut(P,indices,S,scaling_option,optimization_option);
		else
		  min_cut_val=inf; 
		end
		
		min_cut_val_arr(i)=min_cut_val; 
		first_cluster{i}=orig_indices_first_cluster; 
		
	  end
	end
	
	[dummy split_ci]=min(min_cut_val_arr(1:nc)); 
	
	min_cut_val_arr(split_ci)=-1; 
	
	min_cut_val_arr(nc+1)=-1; 
	
	%this splits the cluster. 
	cluster_assignment(first_cluster{split_ci})=nc+1; 
	
	nc=nc+1; 
  end
  
  assignment_kvv=cluster_assignment; 
