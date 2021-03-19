function [f]=plot_hits(map, ds)
% [f]=plot_hits(map, ds)
%
% TODO: <short description>
%
% TODO: <explain vars>
% 
% $Id: plot_hits.m 621 2005-04-27 07:36:15Z dome $
% D. Brugger, 02 February 2005
% viz/plot_hits.m

figure;
% classes
data = ds.data;
[n, m] = size(data);
for k=1:n
  c(k)=str2num(ds.labels{k});
end
ncl = max(c);
idx = ones(ncl,1);
% split ds 
%for k=1:n
%  filtered_data{c(k)}(idx(c(k)),:) = data(k,:);
%  idx(c(k)) = idx(c(k)) + 1;
%end
filtered_data=split_data(ds);
for k=1:ncl
 h{k} = som_hits(map, filtered_data{k});
end
ncl
som_show(map, 'empty', 'Hits');
cmap=hsv(ncl);
%ccmap=flipdim(cmap,1)
ccmap=1-cmap;
som_show_add('hit', cell2mat(h), 'MarkerColor', cmap);
axes('Position', [0.05 1-0.05 0.05 0.05], 'XTick', [], 'YTick', [], ...
    'Visible', 'off');
x = 0.0;
y = 0.0;
for k=1:length(ds.label_names)
  th=text(x, y, ds.label_names{k}, 'Color', ccmap(k,:), ...
	  'BackgroundColor', cmap(k,:));
  e=get(th, 'Extent');
  y = y - 0.1 - e(4);
end
%f = th;
f = gcf;

