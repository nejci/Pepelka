 function assignment=njw_njw_kmeans(S,k)
 % function assignment=njw_njw_kmeans(S,k)
 % Clusters using NJW followed by NJW followed bu K-Means 
  assignment=cluster_spectral_general(S,k,'njw_gen','njw_kmeans'); 
  
  
