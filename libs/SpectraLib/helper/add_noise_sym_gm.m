function Snew=add_noise_sym_gm(S,alpha,seed)
% function Snew=add_noise_sym_gm(S,alpha,seed)
% 
% Adds in the noise to input matrix S so as to "preserve" the signal to
% noise ratio. The noise added is increase the value of S_ij=S_ji by 
% random(0,1)*alpha*sqrt(d_i*d_j)/n.  
  
% use the seed if given else ignore 
  if nargin == 3
	rand('seed',seed); 
  end
  
  n=length(S); 
  d=sum(S,2); % d_i = Sum_j S_ij
  sqrtd=sqrt(d); 
  noise=(alpha/n)*rand(n,n); % uniformly chosen between 0 and alpha/n
  for i=1:n
  	for j=i:n
	  % make the noise propotional to the % signal strength.  
	  noise(i,j) = noise(i,j) *sqrtd(i)*sqrtd(j); 
	  noise(j,i)=noise(i,j); 
	end 
  end
  Snew=S+noise; 
  
  
  
