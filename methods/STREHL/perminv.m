% function ind = perminv(cl)
%
% copyright (c) 1998-2002 by Alexander Strehl

function ind = perminv(cl)
temp = [cl;1:length(cl)];
temp = sortrows(temp',1);
ind = temp(:,2)';
