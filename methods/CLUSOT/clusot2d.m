function [clust,new_sf,new_sl,C,BB,Mx,pMx,BP]=clusot2d(map,td,varargin)
% [clust,sf,sl,C,BB,Mx,pMx,BP]=clusot2d(map,td,varargin)
%
% clusot2d runs the clusot algorithm for arbitrary 2d-topologies.
%
% Inputs:
% map - (som_map_struct) a trained som map
% td - (som_train_struct) training data
%   or (vector) neuron frequencies.
%      If training data is given neuron frequencies
%      are calculated using the given map.
% [argID,   (string) See below. 
%    value] type of value is dependent on argID
% 
% argID can be:
%
% -clusot_2d-------------------------------------------------------------
% cluster_method (string) selector for clustering method, can be
%              either 'clusot_cluster' or 'rec_flood. Default: 'rec_flood'.
%
% -clusot_surf-----------------------------------------------------------
% 'res' - (double) resolution of surface, default 0.05.
% 'type' - (string) computation method 'spline' or 'ellipse'. Default:
%                 'spline'.
% 'ss_flag' - (bool) if set to true, special cases in spline
%             surface computation are handled by inserting mirrored
%             neuron positions. Default: true.
%
% -clusot_cluster--------------------------------------------------------
% 'g' - (double) threshhold parameter in (0,1], default 0.5.
% 'T' - (double) resolution parameter for scanlines, default 256.
% 'flag' - (bool) if set to true, clusters are automerged if
%                 their cluster regions overlap in C. Default: false
%
%
% -rec_flood-------------------------------------------------------------
% 'theta' - (double) waterline threshold, applied w/r to current
%             local maximum, where 0 < theta <= 1. Default: 0.7.
% 'step' - (double) step size for waterline. Default: 0.05.
% 'theta0' - (double) waterline threshold, which is used for
%             preflooding, where 0 <= theta0 < theta. Default: 0.3.
% 'tau' - (double) threshold for controlling the depth of 'valleys'
%             between clusters, where 0 <= tau < 1. Default: 0.
% 'flag' - (bool) indicates whether gradient information is used to
%             adjust cluster borders. Default: false
% 
% Outputs:
% clust - (vector) cluster numbers for map units.
% ------------------------------------------------------------
% sf - (matrix) surface as computed by clusot_surf.m
% sl - (matrix) slope as computed by clusot_surf.m
% ------------------------------------------------------------
% C - (matrix) cluster matrix as computed by clusot_cluster.m
% BB - (matrix) cluster border points (n x 2) matrix
% Mx - (vector) all local maxima
% pMx - (vector) processed local maxima, e.g. those which weren't
%                ignored 
% BP - (matrix) all cluster border points as determined by
%               interpolation method between points in BB (n x 2)
%               matrix
% 
% $Id: clusot2d.m 847 2005-07-25 20:39:02Z dome $
% D. Brugger, 20 March 2005
% algo/clusot2d.m

% Allow caching of surface values?
PRECACHE = 0;

if(nargin == 0)
  test_clusot2d();
  return;
end

% default values
cluster_method='rec_flood';

res=0.05;
type='spline';
ss_flag=true;

T=256;
g=0.5;
flag=false;

theta=0.7; 
step=0.05; 
theta0=0.3;
tau=0;

% parse options
i=1; 
while i<=length(varargin) 
  argok = 1; 
  if ischar(varargin{i}) 
    switch varargin{i} 
      % argument IDs
     case 'cluster_method', i=i+1; cluster_method = varargin{i}; 
     case 'res',       i=i+1; res = varargin{i}; 
     case 'type',      i=i+1; type = varargin{i}; 
     case 'T',         i=i+1; T = varargin{i}; 
     case 'g',         i=i+1; g = varargin{i}; 
     case 'flag',      i=i+1; flag = varargin{i}; 
     case 'ss_flag',   i=i+1; ss_flag = varargin{i}; 
     case 'theta',     i=i+1; theta = varargin{i}; 
     case 'theta0',    i=i+1; theta0 = varargin{i}; 
     case 'step',      i=i+1; step = varargin{i}; 
     case 'tau',       i=i+1; tau = varargin{i}; 
     case 'f',       i=i+1; f = varargin{i}; 
        otherwise
            argok=0;
    end
  else
    argok=0;
  end
  if ~argok, 
    disp(['(clusot2d) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end

%fprintf('clusot2d using:\n');
%fprintf('res = %g\n',res);
%fprintf('type = ''%s''\n', type);
%fprintf('ss_flag = %d\n', ss_flag);
%fprintf('theta = %g\n', theta);
%fprintf('theta0 = %g\n', theta0);
%fprintf('step = %d\n', step);
%fprintf('tau = %g\n', tau);
%fprintf('flag = %d\n',flag);

% use caching
persistent prev_sf;
persistent prev_sl;
persistent prev_sgx;
persistent prev_sgy;
persistent prev_map;
persistent prev_td;
persistent prev_res;
persistent prev_type;
persistent prev_ss_flag;

if(~PRECACHE || isempty(prev_sf) || isempty(prev_sl) || isempty(prev_sgx) || isempty(prev_sgy) ...
   || ~check_matching(map,prev_map) || ~check_matching(td,prev_td) ...
   || ~check_matching(res,prev_res) || ~check_matching(type,prev_type) ...
   || ~check_matching(ss_flag,prev_ss_flag))
  % compute surface + slope
  [new_sf,new_sl,sgx,sgy]=clusot_surf(map,td,res,type,ss_flag);
  % remember values
  prev_sf = new_sf;
  prev_sl = new_sl;
  prev_sgx = sgx;
  prev_sgy = sgy;
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
  sgx = prev_sgx;
  sgy = prev_sgy;
end

% linearly scale slope to [0,1]
new_sl=lin_scale(new_sl,[min(min(new_sl)) max(max(new_sl))],[0 1]);

% compute clusters
if(strcmp(cluster_method,'clusot_cluster'))
  [C,BB,Mx,pMx,BP]=clusot_cluster(new_sf,new_sl,g,T,flag);
elseif(strcmp(cluster_method,'rec_flood'))
  [C,BP,Mx,pMx]=rec_flood(new_sf,sgx,sgy,'theta',theta,'step',step, ...
                          'theta0',theta0,'tau',tau,'flag',flag);
  BB=BP;
else
  error('!!! Unkown cluster method ''%s'' !!!\n',cluster_method)
end

% convert cluster map to clustering
clust=clm2clust(C',map,res,'majority');

