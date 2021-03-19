function L=S_to_L(S) 
% function L=S_to_L(S)
% compute the Laplacian 
% Compute the Stochastic Matrix  P
  d=sum(S,2); 
  n=size(d);
  d=1/sqrt(d); 
  
  for i=1:n
	for j=1:n
	  L(i,j)= S(i,j) *d(i)*d(j); 
	end
  end
  

  
