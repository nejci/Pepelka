function assignment=njw_mcut_ward(S,k)
% function assignment=njw_mcut_ward(S,k)
% Clusters using NJW followed by Multicut and then Ward 

assignment=cluster_spectral_general(S,k,'njw_gen','mcut_ward'); 
  
  
