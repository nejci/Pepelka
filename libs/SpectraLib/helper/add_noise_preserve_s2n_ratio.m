function Snew=add_noise_preserve_s2n_ratio(S,alpha,seed)
% function Snew=add_noise_preserve_s2n_ratio(S,alpha)
% 
% Adds in the noise to input matrix S so as to "preserve" the signal to
% noise ratio. The noise added is increase the value of ROW in S by ratio of alpha
  
  % use the seed if given else ignore 
  if nargin == 3
	rand('seed',seed); 
  end
  
  n=length(S); 
  d=sum(S,2); % d_i = Sum_j S_ij
  noise=(alpha/n)*rand(n,n); % uniformly chosen between 0 and alpha/n
  for i=1:n
    noise(i,:) = noise(i,:) *d(i); % make the noise propotional to the
	                               % signal strength.  
  end 
  Snew=S+noise; 
  
  
  
