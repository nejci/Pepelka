function assignment=sm_ncut(S,k) 
% function assignment=sm_ncut(S,k) 
% Cluster using the ShiMalik algorithm followed bythe normalized cut 

assignment=cluster_shi_r(S,k,'ncut'); 
  
