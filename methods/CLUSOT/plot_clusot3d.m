function []=plot_clusot3d(map,td,varargin)
% []=plot_clusot3d(map,td,varargin)
%
% plot_clusot3d converts a given map for training data td and
% the cluster-matrix C into a scene-graph in open-inventor format
% and stores it in a file.
% This scene-graph is automatically displayed using the viewer application 
%
% [argID,   (string) See below. 
%    value] type of value is dependent on argID
% 
% argID can be:
% C - (double matrix) cluster matrix as returned by
%                     clusotnd. Default: []
% clust - (double vector) clustering as computed by
%                     clm2clust. Default: []
% res - (double) resolution used for surface computation. Default: 0.05
% pMx - (double matrix) processed local maxima as returned by
%                       clusotnd. If given local maxima are
%                       represented by small spheres. Default: []
% cmap - (string) colormap to use, e.g. 'gray'. Default: 'jet'
% fancy - (bool)  whether to use fancy 3d-rendering for hit
%                 numbers. Default: 0 (false)
%
% 
% $Id: plot_clusot3d.m 858 2005-07-25 20:47:10Z dome $
% D. Brugger, 17 June 2005
% viz/plot_clusot3d.m

% default values
res=0.05;
pMx=[];
C=[];
clust=[];
num_n = size(map.codebook,1);
cmap = colormap(jet(num_n));
fancy = 0;
filename = [tempname '.iv'];

% parse options
i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i}, 
      % argument IDs
     case 'res',       i=i+1; res = varargin{i}; 
     case 'pMx',       i=i+1; pMx = varargin{i}; 
     case 'C',         i=i+1; C = varargin{i}; 
     case 'clust',         i=i+1; clust = varargin{i}; 
     case 'cmap',      i=i+1; cmap = cmap_helper(varargin{i}, ...
                                                 num_n);
     case 'fancy',     i=i+1; fancy = varargin{i};
      otherwise argok=0; 
    end
  else
    argok=0;
  end
  if ~argok, 
    disp(['(plot_clusot3d) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end

% increment cluster numbers
C = C + 1;
% create open-inventor file
fid = fopen(filename, 'w');
if(fid == -1)
  error('Could not open file ''%s'' for writing\n', filename);
end
% print header
print_header(fid);

% calculate neuron positions in grid
ndim = size(map.topol.msize,2);
if(ndim <= 2)
  vc = som_vis_coords(map.topol.lattice, map.topol.msize);
  positions = [vc(:,1) abs(vc(:,2)-(max(vc(:,2))))+1 ...
               zeros(size(vc(:,1)))];
elseif(ndim == 3)
  positions = som_unit_coords(map);
%  positions = [vc(:,1) abs(vc(:,2)-(max(vc(:,2))))+1 ...
%               vc(:,3)];
else
  error('!!! Cannot visualize higher dimensional topologies !!!');
end

% determine neuron frequencies
f = calc_f(map,td);
% determine normalized distance matrix between weight vectors
d = calc_d(map,[0,1]);
% for comparison with winnerdist default
cmap = flipdim(cmap,1);

% determine 1-neighborhood in output space
ne1=som_unit_neighs(map.topol);

print_node_begin(fid, 'Separator');
fprintf(fid, 'renderCaching %s\n', 'ON');
print_node_begin(fid, 'DrawStyle');
fprintf(fid, 'lineWidth %d \n', 4);
print_node_end(fid);
print_node_begin(fid, 'Font');
font_sz=4;
fprintf(fid, 'size %d\n', font_sz);
fprintf(fid, 'name "%s\n"', 'Arial');
print_node_end(fid);
if(fancy)
  print_node_begin(fid, 'ProfileCoordinate2');
  fprintf(fid, 'point [ 0 0, 0.05 0.05, 0.25 0.05, 0.3 0 ]\n');
  print_node_end(fid);
  print_node_begin(fid, 'LinearProfile');
  fprintf(fid, 'index [ 0, 1, 2, 3 ]');
  print_node_end(fid);
  print_node_begin(fid, 'Complexity');
  fprintf(fid, 'type %s\n', 'OBJECT_SPACE');
  fprintf(fid, 'value %d\n', 1);
  print_node_end(fid);
end
print_node_begin(fid, 'DirectionalLight');
fprintf(fid,'intensity %g\n', 0.7);
print_node_end(fid);

num_clusters = double(max(max(max(C))));
cmap2 = colormap(jet(num_clusters));

for k=1:num_n
  % draw lines between neuron positions, width proportional
  % to distance
  n_ind = find(ne1(k,:));
  % TODO: put this in aux function draw_sphere_at
  scale_factor=1/res;
  print_node_begin(fid, 'TransformSeparator');
  print_node_begin(fid, 'Translation');
  fprintf(fid, 'translation %d %d %d', scale_factor*positions(k,1), scale_factor*positions(k,2), ...
          scale_factor*positions(k,3));
  print_node_end(fid);
  if(~isempty(clust) && clust(k) > 1)
    print_node_begin(fid, 'BaseColor');
    fprintf(fid, 'rgb %d %d %d\n', cmap2(clust(k),1), cmap2(clust(k),2), ...
            cmap2(clust(k),3));
    print_node_end(fid);
  else
    print_node_begin(fid, 'BaseColor');
    fprintf(fid, 'rgb %d %d %d\n', 1, 1, 1);
    print_node_end(fid);
  end
  
  print_node_begin(fid, 'TransformSeparator');
  print_node_begin(fid, 'Translation');
  fprintf(fid, 'translation %d %d %d', 0, -font_sz/2 ,0);
  print_node_end(fid);
  
  if(fancy)
    print_node_begin(fid, 'Text3');
    fprintf(fid, 'string ["%s"]\n', num2str(f(k)));
    fprintf(fid, 'justification %s\n', 'CENTER');
    fprintf(fid, 'parts %s\n', 'ALL');
    print_node_end(fid);
  else
    print_node_begin(fid, 'AsciiText');
    fprintf(fid, 'justification %s\n', 'CENTER');
    fprintf(fid, 'string %s\n', num2str(f(k)));
    print_node_end(fid);
    print_node_end(fid);
  end

%  print_node_begin(fid, 'Separator');
%  print_node_begin(fid, 'Material')
%  fprintf(fid, 'transparency %g\n', 0.7);
%  fprintf(fid, 'ambientColor %g %g %g\n', 1.0, 0.0, 0.0);
%  fprintf(fid, 'shininess %g\n', 0);
%  print_node_end(fid);
%  print_node_begin(fid, 'BaseColor');
%  fprintf(fid, 'rgb %d %d %d\n', 1, 0, 0);
%  print_node_end(fid);
%  print_node_begin(fid, 'Sphere');
%  fprintf(fid, 'radius %g\n', 2);
%  print_node_end(fid);
%  print_node_end(fid);
  print_node_end(fid); % Transform Separator
  for l=1:length(n_ind)
    % determine line color
    line_color = cmap(max(1,floor(d(k,n_ind(l))*num_n)),:);
    print_node_begin(fid, 'BaseColor');
    fprintf(fid, 'rgb %d %d %d\n', line_color(1), line_color(2), line_color(3));
    print_node_end(fid);
    % draw connections
%    print_node_begin(fid, 'TransformSeparator');
%    height = scale_factor*edist(positions(k,:), ...
%                                positions(n_ind(l),:));
%    % move cylinder to origin
%    print_node_begin(fid, 'Translation');
%    fprintf(fid, 'translation %d %d %d', scale_factor*positions(k,1), ...
%            scale_factor*positions(k,2), ...
%            scale_factor*positions(k,3));
%    print_node_end(fid);

%    % rotate
%    compute_rotation(fid, positions(k,:)+[0 height/2 0], ...
%                     positions(n_ind(l),:)-positions(k,:))
    
%    print_node_begin(fid, 'Cylinder');
%    fprintf(fid, 'radius %d\n', 1);
%    fprintf(fid, 'height %d\n', height);
%    print_node_end(fid);
%    print_node_end(fid);

    print_node_begin(fid, 'Coordinate3');
    % compute new enpoints of connection lines, to
    % avoid overlap with hit numbers
    dx = positions(n_ind(l),1)-positions(k,1);
    dy = positions(n_ind(l),2)-positions(k,2);
    dz = positions(n_ind(l),3)-positions(k,3);
    modf_x = 0; modf_y = 0; modf_z = 0;
    if(dx > 0)
      modf_x = font_sz/2;
    end
    if(dx < 0)
      modf_x = -font_sz/2;
    end
    if(dy > 0)
      modf_y = font_sz/2;
    end
    if(dy < 0)
      modf_y = -font_sz/2;
    end
    if(dz > 0)
      modf_z = font_sz/2;
    end
    if(dz < 0)
      modf_z = -font_sz/2;
    end
    fprintf(fid, 'point [ %d %d %d, %d %d %d ]', ...
            scale_factor*positions(k,1)+modf_x, ...
            scale_factor*positions(k,2)+modf_y, ...
            scale_factor*positions(k,3)+modf_z, ...
            scale_factor*positions(n_ind(l),1)-modf_x, ...
            scale_factor*positions(n_ind(l),2)-modf_y, ...
            scale_factor*positions(n_ind(l),3)-modf_z);
    print_node_end(fid);
    print_node_begin(fid, 'IndexedLineSet');
    fprintf(fid, 'coordIndex [ %d, %d ]\n', 0, 1);
    print_node_end(fid);
  end

end
if(~isempty(C))
  print_point_clouds(fid,C);
end
if(~isempty(pMx))
  for k=1:size(pMx,1)
    print_node_begin(fid, 'BaseColor');
    fprintf(fid, 'rgb %d %d %d\n', 1, 1, 1);
    print_node_end(fid);
    print_node_begin(fid, 'TransformSeparator');
    print_node_begin(fid, 'Translation');
    fprintf(fid, 'translation %d %d %d', pMx(k,1), pMx(k,2), ...
            pMx(k, 3));
    print_node_end(fid);
    print_node_begin(fid, 'Sphere');
    fprintf(fid, 'radius %g\n', 2);
    print_node_end(fid);
    print_node_end(fid);
  end
end
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', -scale_factor, 0, 0);
print_node_end(fid);

print_node_begin(fid, 'TransformSeparator');
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', -scale_factor, 1.5*scale_factor, ...
        -scale_factor);
print_node_end(fid);

print_node_begin(fid, 'TransformSeparator');
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', 1.5, -1, 0);
print_node_end(fid);
print_node_begin(fid, 'BaseColor');
fprintf(fid, 'rgb %d %d %d\n', 1, 1, 1);
print_node_end(fid);
print_node_begin(fid, 'Font');
fprintf(fid, 'size %d\n', 2);
fprintf(fid, 'name "%s\n"', 'Arial');
print_node_end(fid);
print_node_begin(fid, 'AsciiText');
fprintf(fid, 'justification %s\n', 'LEFT');
fprintf(fid, 'string "%s"\n', '0');
print_node_end(fid);
print_node_end(fid);

for k=1:min(num_n,10)
  color = cmap(ceil((k/10)*num_n),:);
  print_node_begin(fid, 'BaseColor');
  fprintf(fid, 'rgb %d %d %d\n', color(1), color(2), color(3));
  print_node_end(fid);
  print_node_begin(fid, 'Cylinder');
  fprintf(fid, 'radius %g\n', 1);
  fprintf(fid, 'height %g\n', 0.5);
  print_node_end(fid);
  print_node_begin(fid, 'Translation');
  fprintf(fid, 'translation %d %d %d', 0, 0.5, 0);
  print_node_end(fid);  
end

print_node_begin(fid, 'TransformSeparator');
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', 1.5, 0, 0);
print_node_end(fid);
print_node_begin(fid, 'BaseColor');
fprintf(fid, 'rgb %d %d %d\n', 1, 1, 1);
print_node_end(fid);
print_node_begin(fid, 'Font');
fprintf(fid, 'size %d\n', 2);
fprintf(fid, 'name "%s\n"', 'Arial');
print_node_end(fid);
print_node_begin(fid, 'AsciiText');
fprintf(fid, 'justification %s\n', 'LEFT');
fprintf(fid, 'string "%s"\n', '1');
print_node_end(fid);
print_node_end(fid);

print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', 0, font_sz, 0);
print_node_end(fid);  
print_node_begin(fid, 'BaseColor');
fprintf(fid, 'rgb %d %d %d\n', 1, 1, 1);
print_node_end(fid);
print_node_begin(fid, 'AsciiText');
fprintf(fid, 'justification %s\n', 'CENTER');
fprintf(fid, 'string ["%s", "%s"]\n', 'Euklidische', 'Distanz');
print_node_end(fid);
print_node_end(fid);

% reset font size
print_node_begin(fid, 'Font');
font_sz=4;
fprintf(fid, 'size %d\n', font_sz);
fprintf(fid, 'name "%s\n"', 'Arial');
print_node_end(fid);

% add small axes object
print_node_begin(fid, 'BaseColor');
fprintf(fid, 'rgb %d %d %d\n', 1, 1, 1);
print_node_end(fid);
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', -scale_factor, -scale_factor, ...
        -scale_factor);
print_node_end(fid);

print_node_begin(fid, 'Sphere');
fprintf(fid, 'radius %g\n', 0.5);
print_node_end(fid);

% x-axis
print_node_begin(fid, 'TransformSeparator');
print_node_begin(fid, 'RotationXYZ');
fprintf(fid, 'angle %g', -pi/2);
fprintf(fid, 'axis %s', 'Z');
print_node_end(fid);
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', 0, scale_factor/2, 0);
print_node_end(fid);
print_node_begin(fid, 'Cylinder');
fprintf(fid, 'radius %g\n', 0.5);
fprintf(fid, 'height %g\n', scale_factor);
print_node_end(fid);
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', 0, scale_factor/2, 0);
print_node_end(fid);
print_node_begin(fid, 'Cone');
fprintf(fid, 'bottomRadius %g', 1); 
fprintf(fid, 'height %g', 1.5);
print_node_end(fid);
print_node_end(fid);
print_node_begin(fid, 'TransformSeparator');
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', scale_factor+font_sz/2, 0, 0);
print_node_end(fid);
print_node_begin(fid, 'AsciiText');
fprintf(fid, 'justification %s\n', 'CENTER');
fprintf(fid, 'string %s\n', 'X');
print_node_end(fid);
print_node_end(fid);

% y-axis
print_node_begin(fid, 'TransformSeparator');
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', 0, scale_factor/2, 0);
print_node_end(fid);
print_node_begin(fid, 'Cylinder');
fprintf(fid, 'radius %g\n', 0.5);
fprintf(fid, 'height %g\n', scale_factor);
print_node_end(fid);
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', 0, scale_factor/2, 0);
print_node_end(fid);
print_node_begin(fid, 'Cone');
fprintf(fid, 'bottomRadius %g', 1); 
fprintf(fid, 'height %g', 1.5);
print_node_end(fid);
print_node_end(fid);
print_node_begin(fid, 'TransformSeparator');
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', 0, scale_factor+font_sz/2, 0);
print_node_end(fid);
print_node_begin(fid, 'AsciiText');
fprintf(fid, 'justification %s\n', 'CENTER');
fprintf(fid, 'string %s\n', 'Y');
print_node_end(fid);
print_node_end(fid);

% z-axis
print_node_begin(fid, 'TransformSeparator');
print_node_begin(fid, 'RotationXYZ');
fprintf(fid, 'angle %g', pi/2);
fprintf(fid, 'axis %s', 'X');
print_node_end(fid);
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', 0, scale_factor/2, 0);
print_node_end(fid);
print_node_begin(fid, 'Cylinder');
fprintf(fid, 'radius %g\n', 0.5);
fprintf(fid, 'height %g\n', scale_factor);
print_node_end(fid);
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', 0, scale_factor/2, 0);
print_node_end(fid);
print_node_begin(fid, 'Cone');
fprintf(fid, 'bottomRadius %g', 1); 
fprintf(fid, 'height %g', 1.5);
print_node_end(fid);
print_node_end(fid);
print_node_begin(fid, 'TransformSeparator');
print_node_begin(fid, 'Translation');
fprintf(fid, 'translation %d %d %d', 0, 0, scale_factor+font_sz/2);
print_node_end(fid);
print_node_begin(fid, 'AsciiText');
fprintf(fid, 'justification %s\n', 'CENTER');
fprintf(fid, 'string %s\n', 'Z');
print_node_end(fid);
print_node_end(fid);

print_node_end(fid); % Separator

% close open-inventor file
fclose(fid);

% call viewer
wd=cd(clusot_root);
arch = computer('arch');
if ~isempty(strfind(arch, 'glnx'))
    system(sprintf('LD_LIBRARY_PATH=/usr/lib:$LD_LIBRARY_PATH %s/viewer_%s %s', clusot_root, arch, filename));
else
    system([clusot_root filesep sprintf('viewer_%s %s', arch, filename)]);
end
cd(wd);
% remove file
delete(filename);

% auxiliary functions
function print_header(fid)
fprintf(fid, '#Inventor V2.1 ascii\n');
function print_node_begin(fid, nodename)
fprintf(fid, '%s {\n', nodename);
function print_node_end(fid)
fprintf(fid, '}\n');
function compute_rotation(fid, v1, v2)
% normalize vectors
v1n = v1 ./ norm(v1);
v2n = v2 ./ norm(v2);
% compute cross product between v1n and v2n,
% e.g. the desired rotation axis
rax = cross(v1n, v2n);
% if all components are zero, we do not perform
% a rotation
rax = round(rax);
if(sum(rax == 0) < 3)
  % compute angle
  angle = acos(dot(v1n, v2n));
  print_node_begin(fid, 'Rotation');
  fprintf(fid, 'rotation %g %g %g %g', rax(1), rax(2), rax(3), mod(angle,pi));
  print_node_end(fid);
end

function print_point_clouds(fid,C)
sz = size(size(C),2);
if(sz > 3)
  error('Cannot plot cluster matrix with dimension > 3');
end

if(sz == 1)
  error('Not supported?\n');
end
if(sz == 2)
end
if(sz == 3)
  num_clusters = double(max(max(max(C))));
  cmap = colormap(jet(num_clusters));
  close;
  fprintf('num_clusters = %d\n', num_clusters);
  for k=2:num_clusters
    % compute coordinates
    [k1,k2,k3] = ind2sub(size(C), find(C==k));
    K = [k1 k2 k3];
    size(K)
    print_node_begin(fid,'Coordinate3');
    fprintf(fid, 'point [\n');
    fprintf(fid,'%d %d %d,\n', K(1:end-1,:)');
    fprintf(fid, '%d %d %d\n', K(end,:)');
    fprintf(fid, ']\n');
    print_node_end(fid);
    % set color
    print_node_begin(fid, 'BaseColor');
    fprintf(fid, 'rgb %g %g %g\n', cmap(k,1), cmap(k,2), cmap(k,3));
    print_node_end(fid);
    print_node_begin(fid, 'PointSet');
    print_node_end(fid);
  end
end


