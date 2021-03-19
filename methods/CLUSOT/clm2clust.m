function [clust]=clm2clust(C,map,res,flag)
% [clust]=clm2clust(C,map,res,flag)
%
% clm2clust converts a cluster-map representation C, as returned by
% clusot_cluster into a clustering clust for the given map with
% resolution res.
%
% Inputs:
% C - (matrix) cluster map
% map - (som_map_struct) trained som map
% res - (double) resolution used by clusot_surf for surface
%         computation
% flag - (string) if set to 'majority' the extent of each neuron is
%         respected when assigning cluster numbers to
%         neurons. Else extent is not respected. Default:
%         'normal' [optional]
%
% Output:
% clust - (vector) clustering, e.g. a vector of size #Neurons x 1
% 
% $Id: clm2clust.m 810 2005-06-21 15:45:19Z dome $
% D. Brugger, 20 March 2005
% util/clm2clust.m

if(nargin == 0)
  test_clm2clust();
  return;
end

if(nargin < 4)
  flag='normal';
end

% determine number of dimensions
ndim = int32(size(size(C),2));
uc=som_unit_coords(map);
% flip in 2d, convention of som-toolbox
if(ndim == 2)
  uc=[uc(:,1) abs(uc(:,2)-(max(uc(:,2))))];
end
ps=mesh_idx(uc,res);


if(strcmp(flag,'majority'))
  msz = map.topol.msize;
  % flip
  tmp = msz(1);
  msz(1) = msz(2);
  msz(2) = tmp;
  extent=(size(C)-1) ./ (msz-1) ./ 2;
  Csz = int32(size(C));
  bound_h = zeros(1,ndim);
  bound_l = zeros(1,ndim);
  for k=1:size(ps,1)
    % compute true extent of neuron, e.g. treat border cases
    for l=1:ndim
      bound_h(l) = ps(k,l)+extent(l);
      if(bound_h(l) > Csz(l))
        bound_h(l) = Csz(l);
      end
      bound_l(l) = ps(k,l)-extent(l);
      if(bound_l(l) < 1)
        bound_l(l) = 1;
      end
      bound_h(l) = int32(bound_h(l));
      bound_l(l) = int32(bound_l(l));
    end

    if(ndim == 2)
      sel=C(bound_l(1):bound_h(1), bound_l(2):bound_h(2));
    end
    if(ndim == 3)
      sel=C(bound_l(1):bound_h(1), bound_l(2):bound_h(2), bound_l(3):bound_h(3));
    end
    if(ndim > 3)
      error('TODO: Not yet implemented for dimensions > 3');
    end
    cn=get_majority(sel);
    clust(k)=cn+1; 
  end
else % do not use extent of neuron
  if(ndim == 2)
    for k=1:size(ps,1)
      clust(k)=C(ps(k,1),ps(k,2))+1;
    end
  end
  if(ndim == 3)
    for k=1:size(ps,1)
      clust(k)=C(ps(k,1),ps(k,2),ps(k,3))+1;
    end
  end
  if(ndim > 3)
    error('TODO: Not yet implemented for dimensions > 3');
  end
end

% make sure cluster numbers form a sequence, e.g. 1 2 1 1 3,
% instead of 1 5 1 1 7
clust=mk_sequence(clust);

function cn=get_majority(sel)
m=max(max(max(sel)));
count=prod(size(sel(find(sel == 0))));
cn=0;
for k=1:m
  new_count=prod(size(sel(find(sel == k))));
  if(new_count > count)
    cn = k;
  end
end

% following aux. functions are from surf_ellip, TODO: factor them out
function b=on_left_border(np,ix)
b = np(1) == ix(1);
function b=on_right_border(np,ix)
b = np(1) == ix(2);
function b=on_lower_border(np,iy)
b = np(2) == iy(1);
function b=on_upper_border(np,iy)
b = np(2) == iy(2);

function test_clm2clust()
map=som_map_struct(4,'msize',[3 4],'rect');
res=0.5;
C=[1 1 1 0 0 0 4; ...
   1 1 0 0 2 0 0; ...
   1 0 0 2 2 2 2; ...
   0 0 0 0 2 2 2; ...
   3 3 0 0 0 2 2];
er=[3 1 1 0 0 1 0 2 0 2 2 4];
er=er+1;
r=clm2clust(C,map,res)
check_equal(r,er,'r','er');

C=[1 1 1 0 0 0 0; ...
   1 1 0 0 2 0 4; ...
   1 0 0 2 2 2 2; ...
   0 0 0 0 2 2 5; ...
   3 3 0 0 0 5 5];
er=[3 1 1 0 0 1 0 2 0 4 2 0];
er=er+1;
r=clm2clust(C,map,res)
check_equal(r,er,'r','er');

fprintf('test_clm2clust succeded\n');
