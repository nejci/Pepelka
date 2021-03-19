function [h]=plot_top(map, td, varargin)
% [h]=plot_top(map, td, varargin)
%
% plot_top visualizes neighborhood structure of a som
% (grid,hex,rand), weight vector distance between neighbor neurons
% and the frequency of each neuron. Works similar to "winnderdist"
% but generalizes its function to arbitrary topologies.
%
% Inputs:
% map - (som_map_struct) a trained som map
% td - (som_train_struct) training data
%   or (vector) neuron frequencies.
% [argID,   (string) See below. 
%    value] type of value is dependent on argID
% 
% argID can be:
% 'cms' - (string) description of colormap, e.g. 'jet',
%         'gray'. Default: 'jet'.
% 'width' - (int) line width used for drawing
%          connections. Default: 2.
% 'ah'   - (axes handle) use axes handle for plotting.
%          Default: Open a new figure.
% Output:
% h - handle to figure
% 
% 
% $Id: plot_top.m 866 2005-07-28 12:09:00Z dome $
% D. Brugger, 01 February 2005
% viz/plot_top.m

% default values
cms='jet';
width=2;
ah=NaN;
cb=true;

% parse options
i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i}, 
      % argument IDs
     case 'cms',      i=i+1; cms = varargin{i}; 
     case 'width',    i=i+1; width = varargin{i}; 
     case 'ah',       i=i+1; ah = varargin{i};
      case 'cb',      i=i+1; cb = varargin{i};
      otherwise, argok=0; 
    end
  else
    argok=0;
  end
  if ~argok, 
    disp(['(plot_top) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end

% check if axes handle was given
if(isnan(ah))
  figure;
  ah = gca;
end

% calculate neuron positions in grid
ndim = size(map.topol.msize,2);
if(ndim <= 2)
  vc = som_vis_coords(map.topol.lattice, map.topol.msize);
  positions = [vc(:,1) abs(vc(:,2)-(max(vc(:,2))))+1 ...
               zeros(size(vc(:,1)))];
elseif(ndim == 3)
  vc = som_unit_coords(map) + 1;
  positions = [vc(:,1) abs(vc(:,2)-(max(vc(:,2))))+1 ...
               vc(:,3)];
else
  error('!!! Cannot visualize higher dimensional topologies !!!');
end

num_n = size(map.codebook,1);
cmap = cmap_helper(cms,num_n);
colormap(cmap);

% compute neuron frequencies, if necessary
if(isstruct(td) && strcmp(td.type, 'som_data'))
  h=calc_f(map,td);
else
  h=td;
end

% determine normalized distance matrix between weight vectors
d = calc_d(map,[0,1]);
% for comparison with winnerdist default
cmap = flipdim(cmap,1);

% text handle vector
th=zeros(1,num_n);
% line handle vector
lh=zeros(1,sum(sum(d > 0)));

set(ah, 'Visible', 'off', 'XLim', [0 max(positions(:,1))], 'YLim', ...
	[0 max(positions(:,2))]);

% determine 1-neighborhood in output space
ne1=som_unit_neighs(map.topol);
pos=1;
for k=1:num_n
  % draw lines between neuron positions, width proportional
  % to distance
  n_ind = find(ne1(k,:));
  for l=1:length(n_ind)
    lh(pos)=line([positions(k,1); positions(n_ind(l),1)], [positions(k,2); ...
		    positions(n_ind(l),2)], [positions(k,3); positions(n_ind(l),3)],'Parent',ah);  
    set(lh(pos), 'Color', cmap(max(1,floor(d(k,n_ind(l))* ...
                                           num_n)),:), 'LineWidth', ...
                 width);
    pos = pos + 1;
  end
end
for k=1:num_n
  th(k) = text(positions(k,1),positions(k,2),positions(k,3),num2str(h(k)),'Parent',ah);
end
set(th, 'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 ...
		    1]);

if cb
% Add a colorbar
hcb = colorbar;
colormap(cmap);
set(hcb, 'OuterPosition', [.025 .31 .07 .34], 'Position', [.03 .348, ...
		    .026 .27]);
% Adjust tick marks
set(hcb, 'YTick', linspace(0,num_n,5), 'YtickLabel', num2str([.0; .25; .5; .75; 1]));
% Add label to colorbar
title(hcb,{'Euclidean','Distance'});
end

h = axes('position',[0.05 0.05 0.01 0.01],'units','normalized');
set(h,'visible','off');
th = text(0,0,map.name);
set(th,'Interpreter','none','FontSize',12);
h = gcf;
set(h, 'Color','w');




