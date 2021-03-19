function [C,BP,Mx,pMx]=rec_flood(sf,sgx,sgy,varargin)
% [C,BP,Mx,pMx]=rec_flood(sf,sgx,sgy,varargin)
%
% rec_flood computes a cluster matrix C and cluster border points
% BP by flooding the surface sf recusively. The flooding is stopped
% when the waterline has reached a value equal to mx*theta, where
% mx is the value of the currently largest local maximum. Redundant
% local maxima are removed by preflooding, controlled by parameter
% theta0. The seperation of cluster regions is controlled by
% parameter tau (see below).
%
% Inputs:
% sf - (matrix) surface function
% sgx - (matrix) gradient dsf/dx
% sgy - (matrix) gradient dsf/dy
% [argID,   (string) See below. 
%    value] type of value is dependent on argID
% 
% argID can be:
% 'theta' - (double) waterline threshold, applied w/r to current
%             local maximum, where 0 < theta <= 1. Default: 0.7.
% 'step' - (double) step size for waterline. Default: 0.05.
% 'theta0' - (double) waterline threshold, which is used for
%             preflooding, where 0 <= theta0 < 1. Default: 0.3.
% 'tau' - (double) threshold for controlling the depth of 'valleys'
%             between clusters, where 0 <= tau < 1. Default: 0.
% 'flag' - (bool) indicates whether gradient information is used to
%             adjust cluster borders. Default: false
%
% Outputs:
% C - (matrix) cluster matrix
% BP - (cell array of matrices) all cluster border points, where
%        BP{k} contains the border points for cluster k given by a 
%        (n x 2) matrix.
% Mx - (vector) all local maxima
% pMx - (vector) processed local maxima, e.g. those which weren't
%                ignored 
% 
% $Id: rec_flood.m 850 2005-07-25 20:41:39Z dome $
% D. Brugger, 18 April 2005
% algo/rec_flood.m

if(nargin == 0)
  test_rec_flood();
  return;
end

mask=imregionalmax(sf);
[r,c]=find(sf.*mask > 0);
% Mx contains x- and y-coordinates
% in the first two columns and the value
% of the local maximum in the third column
Mx=zeros(length(r),3);
for k=1:length(r)
  Mx(k,:)=[c(k) r(k) sf(r(k),c(k))];
end

% sort Mx into descending order (therefore 'flipdim') 
% using value of maximum, e.g. column 3, as key
Mx=flipdim(sortrows(Mx,3),1);

% set default values for parameters
theta=0.7; 
step=0.05; 
theta0=0.3;
tau=0;
flag=false;

% parse options
i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i}, 
      % argument IDs
     case 'theta',       i=i+1; theta = varargin{i}; 
     case 'theta0',      i=i+1; theta0 = varargin{i}; 
     case 'step',        i=i+1; step = varargin{i}; 
     case 'tau',         i=i+1; tau = varargin{i}; 
     case 'flag',        i=i+1; flag = varargin{i}; 
     otherwise argok=0; 
    end
  else
    argok=0;
  end
  if ~argok, 
    disp(['(rec_flood) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end

%fprintf('rec_flood using:\n');
%fprintf('theta = %g\n', theta);
%fprintf('theta0 = %g\n', theta0);
%fprintf('step = %d\n', step);
%fprintf('tau = %g\n', tau);
%fprintf('flag = %d\n',flag);

% do some preflooding
sf(sf <= Mx(1,3)*theta0)=Mx(1,3)*theta0;
% remove local maxima, which were flooded above
idx=2;
for k=2:size(Mx,1)
  if(Mx(k,3) <= Mx(1,3)*theta0)
    break;
  end
  idx = idx + 1;
end
if(idx <= size(Mx,1))
  fprintf('Removed %d local maxima by preflooding\n',size(Mx,1)- ...
          idx+1);
end
Mx=Mx(1:idx-1,:);

% initialize connectivity graph of local maxima
adj=ones(size(Mx,1),size(Mx,1))*Inf;
% set initial region number
rn=1;
% run recursive flooding
C=flood(sf,sgx,sgy,theta,step,Mx,rn,tau,adj,flag);

old_C=C;
C=mk_sequence(old_C);
% note that this is extremely unlikely to happen, but
% check anyway
if(sum(sum(C == old_C)) ~= prod(size(C)))
  tmpfilename=tempname;
  fid=fopen('DAVID_VS_GOLIATH.txt','a');
  fprintf(fid,'theta=%g\n',theta);
  fprintf(fid,'theta0=%g\n',theta0);
  fprintf(fid,'step=%g\n',step);
  fprintf(fid,'tau=%g\n',tau);
  fprintf(fid,'flag=%d\n',flag);
  fprintf(fid,'sf,sgx,sgy have been saved to ''%s''\n',tmpfilename);
  fclose(fid);
  save(tmpfilename,'sf','sgx','sgy');
end

% compute border points for viz with plot_clusot
BP=C2BP(C);

% compute processed maxima for viz with plot_clusot
% TODO: C(Mx(k,2),Mx(k,1)) might be 0 because the active contour
% might have 'jumped' over the local maximum, though this is very
% unlikely to happen, e.g. only for small snakes with <20 points
num_cc=max(max(C));
pMx=zeros(num_cc,3);
for k=1:size(Mx,1)
  val=C(Mx(k,2),Mx(k,1));
  if(val > 0)
    tmp = pMx(val,:);
    if((not(isempty(tmp)) && tmp(1,3) < Mx(k,3)) || ...
       isempty(tmp))
      pMx(val,:)=Mx(k,:);
    end
  end
end

% function for recursive flooding
function [C,rn]=flood(sf,sgx,sgy,theta,step,Mx,rn,tau,adj,flag)
% initialize water line, mx and C
wl=min(min(sf)); mx=Mx(1,3); C=zeros(size(sf));
% treat special case
if(wl == mx || wl >= theta*mx)
  %            ^^^^^^^^^^^^^ This case should not arise but
  % in rare cases does due to numeric issues!? These problems
  % should vanish if variable step truly steps through a list
  % of distince surface heights.
  C = zeros(size(sf));
  C(Mx(1,2),Mx(1,1)) = 1;
  C = rn*C;
%  fprintf('!!! Incrementing region number here !!!\n');
  rn = rn + 1;
%  if(iszero(C))
%    fprintf('!!! iszero(C) rn=%d!!!\n',rn);
%  end
  return;
end

break_flag=1;
while(wl < theta*mx)
  % flood sf upto waterline wl
  sf(sf <= wl) = wl;
  mask = zeros(size(sf));
  mask(sf > wl) = 1;
  L = bwlabel(mask,8);
  
  % update graph edge weights
  pairs=all_ordered_pairs([1:size(adj,1)]);
  for k=1:size(pairs,1)
    if(adj(pairs(k,1),pairs(k,2)) == Inf && ...
       L(Mx(pairs(k,1),2),Mx(pairs(k,1),1)) ~= ...
       L(Mx(pairs(k,2),2),Mx(pairs(k,2),1)))
        weight=min(tau*abs( min(sf(Mx(pairs(k,1),2),Mx(pairs(k,1),1)), ...
                          sf(Mx(pairs(k,2),2),Mx(pairs(k,2),1)))-wl ...
                          ) + wl, ...
                   theta*mx-step);
%       weight=min(tau*abs( min(sf(Mx(pairs(k,1),2),Mx(pairs(k,1),1)), ...
%                           sf(Mx(pairs(k,2),2),Mx(pairs(k,2),1)))-wl ...
%                           ) + wl, ...
%                  theta*min(Mx(pairs(k,1),3),Mx(pairs(k,2),3))-step);
%  fprintf('Setting edge weight (%d,%d)=%g\n',pairs(k,1),pairs(k,2),weight);
      adj(pairs(k,1),pairs(k,2))=weight;
      adj(pairs(k,2),pairs(k,1))=weight;
    end
  end
  % remove flooded edges
  adj(adj <= wl)=0;
  % compute connected components of graph
  [num_cc,cc,cc_idx]=connected_components(adj);
  if(num_cc > 1)
    % partition maxima
    new_Mx=cell(num_cc,1);
    for k=1:num_cc
      sel=Mx(cc_idx{k},:);
      ind=find(sel(:,3) >= wl);
      for l=1:size(sel,1)
        if(isempty(find(ind == l)))
          error(['!!! Detected already flooded local maxima at ' ...
                 '(%d,%d). This should *not* happen !!!'],sel(l,1),sel(l,2));
        end
      end
      new_Mx{k}=sel(ind,:);
    end
    C = zeros(size(C));
%    for k=1:size(new_Mx,1)
%      new_Mx{k}
%    end
    bb=zeros(num_cc,5);
    for k=1:num_cc
      % sort local maxima, descending order
      [tmp,idx]=sortrows(new_Mx{k},3);
      new_Mx{k}=flipdim(tmp,1);
      % adjust adjacency matrix of subgraph
      tmp=cc{k};
      cc{k}=tmp(flipdim(idx,1),:);
      new_Mxk=new_Mx{k};
      % compute union of bounding boxes, as regions already
      % seperated in L might still be connected in the graph
      % depending on the selection of the tau parameter, e.g. for
      % tau > 0
      x_l=Inf; x_h=-Inf; y_l=Inf; y_h=-Inf;
      for l=1:size(new_Mxk,1)
        bbk=bounding_box(L,new_Mxk(l,:));
        x_l=min(bbk(1,1),x_l); x_h=max(bbk(1,2),x_h);
        y_l=min(bbk(2,1),y_l); y_h=max(bbk(2,2),y_h);
      end
      bb(k,:)=[x_l x_h y_l y_h (x_h-x_l)*(y_h-y_l)];
    end
    % sort bounding boxes, in descending order of area size
    [bb,idx]=sortrows(bb,5);
    bb=flipdim(bb,1); idx=flipdim(idx,1);
    % process bounding boxes in order of descending area size
    % to prevent clusters from being lost if they are completely
    % surrounded by another cluster
    for k=1:num_cc
      % mask all other parts of surface
      new_Mxk=new_Mx{idx(k)};
      new_sf=mask_all_others(sf,L,new_Mxk);
%      figure; surf(new_sf);
      bbk=bb(k,:);
      % translate coordinates of local maxima
%      new_Mxk
      new_Mxk(:,1)=new_Mxk(:,1)-bbk(1,1)+1;
      new_Mxk(:,2)=new_Mxk(:,2)-bbk(1,3)+1;
%      new_Mxk
      % recursive flooding
      [Ck,rn]=flood(new_sf([bbk(1,3):bbk(1,4)],[bbk(1,1):bbk(1,2)]), ...
                     sgx([bbk(1,3):bbk(1,4)],[bbk(1,1):bbk(1,2)]), ...
                     sgy([bbk(1,3):bbk(1,4)],[bbk(1,1):bbk(1,2)]), ...
                     theta, step, new_Mxk, rn, tau, cc{idx(k)}, flag);
      %figure; surf(Ck);
      % place result in C at appropriate place
      ind=find(Ck > 0);
      tmp=C([bbk(1,3):bbk(1,4)],[bbk(1,1):bbk(1,2)]);
      tmp(ind)=Ck(ind);
      C([bbk(1,3):bbk(1,4)],[bbk(1,1):bbk(1,2)]) = tmp;
    end 
    break_flag=0; break;
  else
    C = rn*L;
    % use adaptive step size, to prevent local maxima from gettting
    % unintentionally flooded, when step size is too large
    if(size(Mx,1) > 1)
      wl = wl + min(step,abs(Mx(end,3)-wl)-eps);
    else
      wl = wl + step;
    end
  end
end
if(break_flag)
  if(flag)
    % contract boundary using gradient information
    % with the 'active contour'-approach
    % compute vector field
    [X,Y]=meshgrid(1:size(sf,2),1:size(sf,1));
    v_x=zeros(size(X)); v_y=zeros(size(Y));
    v_norm=ones(size(X)).*Inf;
    hf=ones(size(X));
    for k=1:size(Mx,1)
      mx_x = Mx(k,1); mx_y = Mx(k,2);
      v_x_next = mx_x-X; v_y_next = mx_y-Y;
      v_norm_next = sqrt(v_x_next.^2+v_y_next.^2);
      ind=find(v_norm_next <= v_norm);
      % disturb to prevent div by zero
      v_norm_next(v_norm_next == 0)=eps;
      v_norm(ind) = v_norm_next(ind);
      v_x(ind) = v_x_next(ind) ./ v_norm_next(ind); 
      v_y(ind) = v_y_next(ind) ./ v_norm_next(ind);
      hf(ind)=abs(sf(ind)-Mx(k,3));
    end
    % figure; quiver(v_x,v_y);
    % compute dot product
    dp = sgx .* v_x + sgy .* v_y;
    % normalize
    dp = lin_scale(dp,[min(min(dp)) max(max(dp))],[0 1]);
    % compute normalized height factor
    hf = lin_scale(hf,[min(min(hf)) max(max(hf))],[0 1]);
    % compute external energy
    fx = (1-dp).*v_x.*hf;
    fy = (1-dp).*v_y.*hf;
    figure; quiver(fx,fy);
    BP=find_contour_helper(C,Mx);
    
    % TODO: might want to lift active contour parameters
    % deform active contour
    [x,y]=snakeinterp1(BP(:,1),BP(:,2),0.5);
    [x,y]=snakedeform(x',y',0.05,0.09,0.05,0.02,fx,fy,20, ...
                      [1, size(sf,2)], [1, size(sf,1)]);
    BP=[round(x) round(y)];
    
    % compute new C
    C=zeros(size(C));
    BP=unique(BP,'rows');
    % interpolate to ensure closed borders
    BP=lin_ipol(BP);
    C=vector_edge_fill(C,BP,0,1);
    C=rn*C;
  end

  % finally, increment region number
  rn = rn + 1;
%  if(iszero(C))
%    fprintf('!!! iszero(C) rn=%d!!!\n',rn);
%  end
%  if(max(max(L))>1)
%    fprintf('!!!! SHOULD NOT BE =%d !!!!\n',max(max(L)));
%  end
%  figure; surf(C);
%  fprintf('!!! Incremented region number to %d !!!\n',rn);
end

function b=iszero(M)
[n,m]=size(M);
b = sum(sum(zeros(n,m) == M)) == n*m;

% mask all areas which do not belong to local maxima Mx
function new_sf=mask_all_others(sf,L,Mx)
new_sf=zeros(size(sf));
for k=1:size(Mx,1)
  idx=L(Mx(k,2),Mx(k,1));
  ind=find(L == idx);
  new_sf(ind)=sf(ind);
end

% auxiliary function used by bounding_box and flood
function K=find_contour_helper(L,Mx)
mx=Mx(1,1:2);
% determine point p with 4-neighbor not in region, w/r to mx
p=[];
if(mx(1) == 1)
  p=mx;
else
  for k=mx(1)-1:-1:1
    if(L(mx(2),k) == 0)
      p=[k+1 mx(2)];
      break;
    end
  end
  % might be the case if the region extends all the way to the left of
  % the plane
  if(isempty(p))
    p=[1 mx(2)];
  end
end
% determine contour
K=find_contour(L,p,L(mx(2),mx(1)));

function bb=bounding_box(L,Mx)
K=find_contour_helper(L,Mx);
% compute bounding box
bb=[min(K(:,1)) max(K(:,1)); ...
    min(K(:,2)) max(K(:,2))];

function test_rec_flood()
map_path='wine_train_lininit_batch_11x6-rect_gaussian_16_60_inv.mat';
load_wine
dr=clusot_root();
wd=cd(sprintf('%s/results/som/wine',dr));
load(map_path);
[sf,sl,sgx,sgy]=clusot_surf(sM,wine_train,0.05,'ellipse');
[C,BP,Mx,pMx]=rec_flood(sf,sgx,sgy,'flag',true);
plot_clusot(sf,sl,C,BP,pMx);
cd(wd);
