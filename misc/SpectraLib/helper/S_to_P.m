function P=S_to_P(S) 
% function P=S_to_P(S) 
% Compute the Stochastic Matrix  P
  P=S;
  D=sum(S,2); 
  n=size(D); 
  for i=1:n
	P(i,:) =P(i,:)/D(i); 
  end
  

  
