function BP=C2BP(C)
% BP=C2BP(C)
%
% compute boundary points of regions in cluster matrix C
%
% Input:
% C - (matrix) cluster matrix.
%
% Output:
% BP - (cell array) boundary points of cluster regions.
%
% $Id: C2BP.m 940 2005-11-04 17:57:33Z dome $
% D. Brugger, 28 July 2005
% util/C2BP.m

BP={};
num_cc=max(max(C));
for k=1:num_cc
  [ys,xs]=find(C == k);
  % compute point p with 4 neighbor not in current region
  [val,idx]=min(xs);
  p=[xs(idx),ys(idx)];
  % compute contour of region
  BP{k}=find_contour(C,p,k);
end
