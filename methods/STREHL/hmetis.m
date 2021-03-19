% function labels=hmetis(x,k,w)
%
% copyright (c) 1998-2002 by Alexander Strehl

function labels=hmetis(x,k,w) 

if ~exist('w','var'),
  filename = wgraph(x,[],2);
else
  filename = wgraph(x,w,3);
end; 
labels = sgraph(k,filename); 
delete(filename);
