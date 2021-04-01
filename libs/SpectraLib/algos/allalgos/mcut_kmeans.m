function assignment=mcut_kmeans(S,k)
% function assignment=mcut_kmeans(S,k)
% Clusters using multicut followed by K-means

  assignment=cluster_spectral_general(S,k,'mcut_all','kmeans'); 
  
  
