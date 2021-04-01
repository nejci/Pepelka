function sig_str=sigma_str(sigma) 
% function sig_str=sigma_str(sigma)
% String representation of given number. Used to make the filename. 
  sig_str=sprintf('%1.2f',sigma); 
%   if length(sig_str)==2
% 	sig_str=strcat('0',sig_str); 
%   end
  
  % added this to make filenames smaller 
  if sigma<0 
	sig_str=''; 
  end
  
  
