function [h]=plot_hit_num(sM,td)
% [h]=plot_hit_num(sM,td)
%
% plot_hit_num creates a hit plot for each class
% seperately. With this kind of plot it is easier
% to judge the quality of a map w/r to class seperation.
%
% Inputs:
% sM - (som_map_struct) a trained som map
% td - (som_train_struct) training data
%
% Ouput:
% h - handle to current figure
% 
% $Id: plot_hit_num.m 864 2005-07-28 11:53:56Z dome $
% D. Brugger, 22 May 2005
% viz/plot_hit_num.m

h=figure;
% classes
data = td.data;
[n, m] = size(data);
for k=1:n
  c(k)=str2num(td.labels{k});
end
ncl = max(c);

% split data according to class labels
filtered_data=split_data(td);

% compute subplot size
n=ceil(sqrt(ncl)); m=ceil(ncl/n);

% use plot_top with optional handle argument to
% generate a hit plot for each class
for k=1:ncl
  h=subplot(n,m,k);
  title(h,sprintf('Class #%d',k));
  plot_top(sM, filtered_data{k}, 'cms', 'gray', 'ah', h, 'cb', false);
end


