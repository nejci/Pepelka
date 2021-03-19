function assignment=njw_kmeans(S,k)
% function assignment=njw_kmeans(S,k)
% Cluster using the NJW followed by the K-means 

assignment=cluster_spectral_general(S,k,'njw_gen','kmeans'); 
  
  
