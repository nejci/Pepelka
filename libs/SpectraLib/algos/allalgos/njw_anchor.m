function assignment=njw_anchor(S,k)
% function assignment=njw_anchor(S,k)
% Cluster using the NJW algorithm followed by the Anchor 
  assignment=cluster_spectral_general(S,k,'njw_gen','anchor'); 
  
  
