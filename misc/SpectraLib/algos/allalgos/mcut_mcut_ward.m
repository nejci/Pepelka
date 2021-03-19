function assignment=mcut_mcut_ward(S,k)
% function assignment=mcut_mcut_ward(S,k)
% Cluster using multicut (two stages) followed by Ward

assignment=cluster_spectral_general(S,k,'mcut_all','mcut_ward'); 
  
  
