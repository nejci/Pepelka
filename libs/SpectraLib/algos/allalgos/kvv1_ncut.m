function assignment=kvv1_ncut(S,k) 
% function assignment=kvv1_ncut(S,k) 
% Clusters using KVV minimizing normalized cut . 
assignment=cluster_kvv(S,k,'scale_row','ncut'); 
  
