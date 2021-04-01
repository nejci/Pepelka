function [assignment_linkage,param_string]=cluster_linkage(S,k,link_method)
  
  % Clusters the input similarity based on the single linkage as
  % implemented in the statistics toolbox in matlab. 
  %  
  % INPUT:
  % S = similarity matrix (called the Affinity Matrix in the paper)
  % k = number of clusters
  %
  % OUTPUT
  % assignment_linkage( 1,n ) vector of integers 1:k assigning the points to clusters


  n=length(S); 
  % convert S to an n(n-1)/2 vector as is the output pclust and is the
  % input desired by linkage
  
  sim_vector=zeros(1,n*(n-1)/2) ;
  index=0; 
  for i=1:(n-1)
        j=(i+1):n;
        num_j=length(j);
        sim_vector(1,index+(1:num_j))=1./(S(i,j)+1e-10);
        index=index+num_j;
     
  end
  
  % compute the single linkage and produce the hierarchical tree
  htree=linkage(sim_vector,link_method); 
  
  % do the clustering 
  assignment_linkage=cluster(htree,k); 
  if nargout >1 
	param_string='';
  end
