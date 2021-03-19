function [clust,new_sf,new_sl,C,Mx,pMx,BP]=clusotnd(map,td,varargin)
% [clust,sf,sl,C,Mx,pMx,BP]=clusotnd(map,td,varargin)
%
% clusotnd runs the clusot algorithm for arbitrary 2d-topologies
% and grid topologies with n dimensions.
%
% Inputs:
% map - (som_map_struct) a trained som map
% td - (som_train_struct) training data
%   or (vector) neuron frequencies.
% [argID,   (string) See below. 
%    value] type of value is dependent on argID
% 
% argID can be:
%
% -clusot_surf-----------------------------------------------------------
% 'res' - (double) resolution of surface, default 0.05.
% 'type' - (string) computation method 'ellipse' or 'spline'. Default:
% 'ss_flag' - (bool) if set to true, special cases in spline
%             surface computation are handled by inserting mirrored
%             neuron positions. Default: true.
%
% -rec_flood-------------------------------------------------------------
% 'theta' - (double) waterline threshold, applied w/r to current
%             local maximum, where 0 < theta <= 1. Default: 0.7.
% 'theta0' - (double) waterline threshold, which is used for
%             preflooding, where 0 <= theta0 < theta. Default: 0.3.
% 'tau' - (double) threshold for controlling the depth of 'valleys'
%             between clusters, where 0 <= tau < 1. Default: 0.
%
% Outputs:
% clust - (vector) cluster numbers for map units.
% ------------------------------------------------------------
% sf - (matrix) surface as computed by clusot_surf.m
% sl - (matrix) slope as computed by clusot_surf.m
% ------------------------------------------------------------
% C - (matrix) cluster matrix as computed by clusot_cluster.m
% Mx - (vector) all local maxima
% pMx - (vector) processed local maxima, e.g. those which weren't
%                ignored 
% BP - (cell array) cluster border points. Only computed in 2 dimesnions.
%
% $Id: clusotnd.m 862 2005-07-25 20:50:35Z dome $
% D. Brugger, 01 June 2005
% algo/clusotnd.m

% default values
res=0.05;
type='spline';
ss_flag=true;

theta=0.7; 
theta0=0.3;
tau=0;
flag=0;

% parse options
i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i}, 
      % argument IDs
     case 'res',       i=i+1; res = varargin{i}; 
     case 'type',      i=i+1; type = varargin{i}; 
     case 'ss_flag',   i=i+1; type = varargin{i}; 
     case 'flag',      i=i+1; flag = double(varargin{i}); 
     case 'theta',     i=i+1; theta = varargin{i}; 
     case 'theta0',    i=i+1; theta0 = varargin{i}; 
     case 'tau',       i=i+1; tau = varargin{i}; 
     otherwise
         argok=0;
    end
  else
    argok=0;
  end
  if ~argok, 
    disp(['(clusotnd) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end

ndim = size(map.topol.msize,2);
% adjust surface options
if(ndim > 2)
  type = 'ellipse';
end

% use caching
persistent prev_sf;
persistent prev_sl;
persistent prev_map;
persistent prev_td;
persistent prev_res;
persistent prev_type;
persistent prev_ss_flag;

if(isempty(prev_sf) || isempty(prev_sl) ...
   || ~check_matching(map,prev_map) || ~check_matching(td,prev_td) ...
   || ~check_matching(res,prev_res) || ~check_matching(type,prev_type) ...
   || ~check_matching(ss_flag,prev_ss_flag))
  
  % compute surface + slope
  if(strcmp(type,'spline') || ndim == 2 || ...
     strcmp(map.topol.lattice,'hexa'))
    % use matlab implementation
    tic; [new_sf,new_sl]=clusot_surf(map,td,res,type,ss_flag); t=toc;
    fprintf('Elapsed time for clusot_surf is %g seconds.\n', t);	
  else
    % use c implementation
    
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
    if(ndim < 3)
      % mirror coordinate in y-axis direction, so surface layout can be
      % directly compared to output of plot_top etc.
      uc=[uc(:,1) abs(uc(:,2)-(max(uc(:,2))))];
    end

    % adjust msize to reflect flipping of x- and y-coordinates, which is
    % done by som_unit_coords.
    msize = map.topol.msize;
    msize = [msize(2) msize(1) msize(3:end)];
    % compute surface
    tic; [new_sf,new_sl]=c_surf_ellip(map.codebook, hn, uc, msize, res); t=toc;
    fprintf('Elapsed time for c_surf_ellip is %g seconds.\n', t);	

    if(ndim < 3)
      % transpose, so surface layout can be compared to output of plot_top
      new_sf=new_sf';
      new_sl=new_sl';
    end
  end
  
  % remember values
  prev_sf = new_sf;
  prev_sl = new_sl;
  prev_map = map;
  prev_td = td;
  prev_res = res;
  prev_type = type;
  prev_ss_flag = ss_flag;
  
else
  % else use cached values
  fprintf('Using cached surface values...\n');
  new_sf = prev_sf;
  new_sl = prev_sl;
end

% linearly scale slope to [0,1]
new_sl=lin_scale(new_sl,[min(min(new_sl)) max(max(new_sl))],[0 1]);

% compute clusters
tic; 
[C,Mx,pMx]=c_rec_flood(new_sf, 'theta', theta, 'theta0', theta0, 'tau', tau, 'flag', flag); 
t=toc;
fprintf('Elapsed time for c_rec_flood is %g seconds.\n', t);

if(ndim < 3)
  pMx = [pMx(:,2) pMx(:,1) pMx(:,3)];
  Mx = [Mx(:,2) Mx(:,1) Mx(:,3)];
end

C = double(C); C = C - 1;
% convert cluster map to clustering
clust=clm2clust(C,map,res,'majority');

% compute boundary points, if possible
if(ndim < 3)
  BP=C2BP(C);
else
  BP={};
end


