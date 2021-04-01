function assignment=mcut_anchor(S,k)
% function assignment=mcut_anchor(S,k)
% Clusters using multicut followed bu anchor algorithm 
  assignment=cluster_spectral_general(S,k,'mcut_all','anchor'); 
  
  
