function [sf,sl,sgx,sgy]=clusot_surf(map,td,res,type,bflag)
% [sf,sl,sgx,sgy]=clusot_surf(map,td,res,type,bflag)
%
% clusot_surf computs the surface function sf and slope sl
% for a trained som map map for training data td. res and
% type determine resolution of surface and the computation
% method used.
%
% Inputs:
% map - (som_map_struct) trained som map
% td - (som_train_struct) training data
%   or (vector) neuron frequencies.
%      If training data is given neuron frequencies
%      are calculated using the given map.
% res - (double) resolution of surface
% type - (string) computation method 'spline' or 'ellipse'
%         Default for grid topology is 'ellipse' for other
%         topologies default is 'spline'
% bflag - (bool) if set to true the handling of special cases
%                 uses mirrored neuron positions. Default:
%                 true. [optional]
%
% Outputs:
% sf - (matrix) surface function
% sl - (matrix) slope
% 
% $Id: clusot_surf.m 849 2005-07-25 20:40:13Z dome $
% D. Brugger, 14 March 2005
% algo/clusot_surf.m

if(nargin == 0)
  test_clusot_surf();
  return;
end

if(nargin < 4)
  type='ellipse';
end

if(nargin < 5)
  bflag=true;
end

% determine dimension of euclidian plane ...
msz=map.topol.msize;
%ix=[1 msz(2)]; iy=[1 msz(1)];
ix=[0 msz(2)-1]; iy=[0 msz(1)-1];
% and total number of neurons
num_n=msz(1)*msz(2);

% compute normalized distances [0,0.99] 
% between codebook vectors
d=calc_d(map,[0,0.99]);

% compute neuron frequencies, if necessary
if(isstruct(td) && strcmp(td.type, 'som_data'))
  h=calc_f(map,td);
else
  h=td;
end
% normalize frequencies
hn=h./sqrt(2*pi);
  

% compute neuron positions
uc=som_unit_coords(map);
% mirror coordinate in y-axis direction, so surface layout can be
% directly compared to output of plot_top etc.
uc=[uc(:,1) abs(uc(:,2)-(max(uc(:,2))))];

% neighborhood matrix
ne=som_unit_neighs(map.topol);
xdim=round(ix(2)/res-ix(1)/res+1); ydim=round(iy(2)/res-iy(1)/res+1);
sf=zeros(ydim,xdim);
sgx=zeros(ydim,xdim); sgy=zeros(ydim,xdim);

if(strcmp(map.topol.lattice,'hexa') || ...
   strcmp(type, 'spline'))
   % for all neurons in map
   for k=1:num_n
     % position of current neuron N_ij
     np=uc(k,:);
     % compute indices of neighbor units
     ind=find(ne(k,:));
     p=zeros(size(ind,2),2);
     % compute positions of control points
     for l=1:size(ind,2)
       % p = (m-N_ij)*(1-d(N_ij,m))+N_ij
%       uc(ind(l),:)
%       np
%       d(k,ind(l))
       p(l,:)=(uc(ind(l),:)-np)*(1-d(k,ind(l)))+np;
     end
     % compute surface and partial derivatives
     [f,dfdx,dfdy]=surf_spline(np,p,hn(k),ix,iy,res,bflag);
     % compute sum of surface functions and gradient
%      if(k >= 16 && k <= 20)
%      figure; surf(f); xyz_label();
%      caxis([0 12]);
%      view(-20,50);
%      end
     sf = sf + f;
     sgx = sgx + dfdx;
     sgy = sgy + dfdy;
%     figure; surf(f); xyz_label(); title(sprintf('step %d',k));
%     dump=input(sprintf('Completed step %d. Type return to continue...', ...
%                        k));
%     if(k==36)
%       return;
%     end
   end
   % compute slope
   sl = sqrt(sgx.^2+sgy.^2);
else
   for k=1:num_n
     % position of current neuron N_ij
     np=uc(k,:);
     % compute surface and partial derivatives
     [f,dfdx,dfdy]=surf_ellip(np,d,hn(k),ix,iy,res);
     % compute sum of surface functions and gradients
     sf = sf + f;
     sgx = sgx + dfdx;
     sgy = sgy + dfdy;
%     if(k >= 27 && k <= 36)
%       figure; surf(f); xyz_label();
%     end
%     dump=input(sprintf('Completed step %d. Type return to continue...', ...
%                        k));
%     if(k==36)
%       return;
%     end
   end
   sl = sqrt(sgx.^2+sgy.^2);
end


function test_clusot_surf()
test_dataset='wine';
test_map='train_lininit_batch_4x4-rect_gaussian_4_16_inv';
res=0.05;
type='spline';
[s,clusot_root]=system('echo $DIPLOM_ROOT');
% strip trailing newline
clusot_root=clusot_root(1:end-1);
test_map_path=sprintf('%s/results/som/%s/%s_%s.mat', ...
                      clusot_root, test_dataset, test_dataset, test_map)
load(test_map_path);
eval(sprintf('load_%s',test_dataset));
[sf,sl]=clusot_surf(sM,eval(sprintf('%s_train',test_dataset)),res,type);
figure; surf(sf); xyz_label(); title(sprintf('Surface %s', type));
figure; surf(sl); xyz_label(); title(sprintf('Slope %s', type));
