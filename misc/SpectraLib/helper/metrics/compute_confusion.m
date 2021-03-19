function [n, k_true, k_all,confusion]=compute_confusion(true_assign,assign)
% function [n, k_true, k_all,confusion]=compute_confusion(true_assign,assign)
% Compute the relative confusion matrix for the comparitive assignment. 

  n=length(find (true_assign > 0) );
  k_true=max(true_assign);
  k_all=max(assign);
  k_all = max(k_all, k_true); 
  
  if size(true_assign,2)==1 % single column mode
	true_assign=true_assign';
  end
  if size(assign,2)==1 % single column mode
	assign=assign';
  end
  
  confusion = zeros( k_true, k_all );     % build confusion matrix
  for i1 = 1:k_true;
	for i2 = 1:k_all;
	  confusion( i1, i2 ) = length( find( (true_assign == i1) & (assign == i2 )));
	end;
  end;

%  confusion = confusion/n;
