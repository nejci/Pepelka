function assignment=njw_mcut_kmeans(S,k)
% function assignment=njw_mcut_kmeans(S,k)
% Clusters using multicut followed by kmeans. 
  assignment=cluster_spectral_general(S,k,'njw_gen','mcut_kmeans'); 
  
  
