function assignment=kvv1_cond(S,k) 
% function assignment=kvv1_cond(S,k)
% Clusters using KVV1 with conductance 

  assignment=cluster_kvv(S,k,'scale_row','conductance'); 
  
