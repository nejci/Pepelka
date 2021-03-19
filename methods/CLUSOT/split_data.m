function [r]=split_data(ds)
% [r]=split_data(ds)
%
% split_data splits matrix ds.data into submatrices
% which only contain patterns which belong to the same
% class, e.g. two classes resuls in to submatrices.
%
% Input:
% ds - (som_data_struct) data to be split
% 
% Output:
% r  - (cell_array) r{k} contains k-th som_data_struct
%
% $Id: split_data.m 390 2005-02-07 21:56:43Z dome $
% D. Brugger, 06 February 2005
% util/split_data.m

%if(nargin == 0)
%  split_data_test();
%  return;
%end

if(isempty(ds.data))
  r = {};
  return;
end

% classes
data = ds.data;
[n, m] = size(data);
for k=1:n
  c(k)=str2num(ds.labels{k});
end
ncl = max(c);
idx = ones(ncl,1);
% init r
for k=1:ncl
  r{k}.type = 'som_data';
  r{k}.data = [];
  r{k}.labels = [];
  r{k}.comp_names = ds.comp_names;
  r{k}.name = ds.name;
  r{k}.comp_norm = ds.comp_norm;
  r{k}.label_names = ds.label_names;
end
% split ds 
for k=1:n
  r{c(k)}.data(idx(c(k)),:) = data(k,:);
  r{c(k)}.labels{idx(c(k)),1} = ds.labels{k};
  idx(c(k)) = idx(c(k)) + 1;
end

function split_data_test()
ds.data = [1 2 3; 4 5 6; 7 8 9; 10 11 12];
ds.labels = {'1';'2';'1';'3'};
ds.name = 'test.dt';
ds.comp_names = {};
ds.comp_norm = {};
ds.label_names = {'c1'; 'c2'; 'c3'};
r=split_data(ds);
if(length(r) ~= 3)
  error('!!! length(r) = %d, but should be %d !!!', length(r), 3);
end
[n,m]=size(r{1}.data);
if(n ~= 2)
  error('!!! length(r{1}.data) = %d, but should be %d !!!', n, 2);
end
[n,m]=size(r{2}.data);
if(n ~= 1)
  error('!!! length(r{2}.data) = %d, but should be %d !!!', n, 1);
end
[n,m]=size(r{2}.data);
if(n ~= 1)
  error('!!! length(r{3}.data) = %d, but should be %d !!!', n, 1);
end
%r{1}.data
%r{2}.data
%r{3}.data
if( sum(sum(r{1}.data == [1 2 3; 7 8 9])) ~= 6 || ...
  sum(sum(r{2}.data == [4 5 6])) ~= 3 || ...
  sum(sum(r{3}.data == [10 11 12])) ~= 3)
  error('!!! Did not get expected contents !!!', length(r{3}), 1);
end  
