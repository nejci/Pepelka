function ari=adjusted_rand_index(true_assign,assign)
% function ari=adjusted_rand_index(true_assign,assign)
%
% Calculates the rand index given the true_assignment and the
% assignment. 
  [n, k_true, k_all,confusion]=compute_confusion(true_assign,assign);
  ari=adjust_rand_index_c(confusion); 
  
	
  
  