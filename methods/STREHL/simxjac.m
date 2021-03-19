% function s = simxjac(a,b)
%
% DESCRIPTION
%   computes extended Jaccard similarity between row objects in matrices a and b
%
% Copyright (c) 1998-2002 by Alexander Strehl

function s = simxjac(a,b)

if ~exist('b','var')
  b = a;
end;

n = size(a,1);
m = size(b,1);
d = size(a,2);
if (d~=size(b,2))
  disp('simxjac: data dimensions do not match');
  return;
end;

temp = double(a) * double(b'); %spremenil Nejc - prej so bile tu logical vrednosti.
asquared = sum((a.^2),2); 
bsquared = sum((b.^2),2); 
s = temp ./ ((asquared * ones(1,m)) + (ones(n,1) * bsquared') - temp);


