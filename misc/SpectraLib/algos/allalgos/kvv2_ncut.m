function assignment=kvv2_ncut(S,k) 
% function assignment=kvv2_ncut(S,k) 
% Clusters using KVV2 algorithm minimizing the normalized cut. 
  assignment=cluster_kvv(S,k,'add_diagonal','ncut'); 
  
