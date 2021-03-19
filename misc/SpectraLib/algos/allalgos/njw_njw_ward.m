function assignment=njw_njw_ward(S,k)
% function assignment=njw_njw_ward(S,k)
% Cluster using NJW followed by NJW followed by Ward 

assignment=cluster_spectral_general(S,k,'njw_gen','njw_ward'); 
  
  
