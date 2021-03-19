function [f]=plot_pca(td,dim,labels)
% [f]=plot_pca(td,dim,labels)
%
% plot_pca creates a pca-plot of dimension dim <= 3
% for the dataset split classwise, e.g. varargin{1}=patterns of
% class1, etc.
%
% Input:
% td  - (som_struct) training data
% dim - (int) dimension, usually 2 or 3
% labels - (vector) training data labels, overrides td.labels [optional]
%
% Output:
% f - handle to current figure
% 
% $Id: plot_pca.m 859 2005-07-25 20:47:58Z dome $
% D. Brugger, 07 February 2005
% viz/plot_pca.m


if(dim > 3)
  error('!!! Doesn''t make sense trying to plot %d dimensions', ...
	dim);
end

n=size(td.data,1);
if(nargin < 3)
  % determine colors
  for k=1:n
    labels(k)=str2num(td.labels{k});
  end
end

%if(max(labels) < 2)
%  error('!!! Problem with less than 2 classes? !!!');
%end

% mark unassigned data
max(labels)
labels(labels == 0) = max(labels)+1;
cmap = hsv(max(labels));
cmap(max(labels),:) = [0 0 0];
colD = cmap(labels,:);
%for k=1:n
%  colD(k,:)=cmap(labels(k),:);
%end
markermap = {'o','+','x','*','v','^','<','>','h','s','d','p','.'};
%markermap = {'o','+','x','.','v','^','<','>','h','s','d','p','.'};
markerD = {markermap{labels}}';

% pca projection
Pd=pcaproj(td.data,dim);

figure;
som_grid('rect',[n, 1],'Line','none','Coord', Pd,'markercolor',colD, ...
         'marker', markerD)
title(sprintf('%dD-PCA of %s data', dim, td.name));
xlabel('pc1');
if(dim >= 2)
  ylabel('pc2');
end
if(dim == 3)
  zlabel('pc3');
end

% add color legend
axes('Position', [0.05 1-0.05 0.05 0.05], 'XTick', [], 'YTick', [], ...
    'Visible', 'off');
x = 0.0;
y = 0.0;
for k=1:max(labels)
  th=text(x, y, 'hallodo', 'Color', cmap(k,:), ...
	  'BackgroundColor', cmap(k,:));
  e=get(th, 'Extent');
  y = y - 0.1 - e(4);
end

f=gcf;




