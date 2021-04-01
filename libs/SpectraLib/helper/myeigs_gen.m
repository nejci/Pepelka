function [v,d,return_indices]=myeigs_gen(A,B,k) 
% function [v,d]=myeigs(S,k) 
% Generalized EV of Ax=lambda*B*x. Uses eigs. 
% Note: d returned is a vector not a diag matrix 
  opts.disp=0; 
  [v d]=eigs(A,B,k,'la',opts); 
  d=diag(d); 
  [dval indices]=sort(-d); 
  dval=-dval ; 
  v=v(:,indices); 
  
  if nargin==3
	return_indices=indices; 
  end
  
