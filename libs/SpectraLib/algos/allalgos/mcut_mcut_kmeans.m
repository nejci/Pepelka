function assignment=mcut_mcut_kmeans(S,k)
% function assignment=mcut_mcut_kmeans(S,k)
% Cluster using multicut (2 stages) followed by K-means 

assignment=cluster_spectral_general(S,k,'mcut_all','mcut_kmeans'); 
  
  
