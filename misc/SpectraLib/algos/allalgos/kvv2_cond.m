function assignment=kvv2_cond(S,k) 
% function assignment=kvv2_cond(S,k) 
% Clusters using the second kvv2 algorithms minimizing conductance. 

  assignment=cluster_kvv(S,k,'add_diagonal','conductance'); 
  
