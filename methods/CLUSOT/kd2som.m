function [sM,f]=kd2som(map_file,hit_file)
% [sM,f]=kd2som(map_file,hit_file)
%
% Utility to build up som_map_struct from kd2 data
% given in map and hit files.
%
% Input:
% map_file - (string) file with codebook vectors
% hit_file - (string) file with map hits
%
% Output:
% sM - (som_map_struct) 
% f - (vector) neuron frequencies
%
% $Id: kd2som.m 846 2005-07-25 20:37:48Z dome $
% D. Brugger, 25 July 2005
% util/kd2som.m

if(nargin == 0)
  test_kd2som();
  return;
end

cbv = load(map_file);
hits = load(hit_file);
[n,m] = size(hits);

% number of neurons
num_n = n*m;
% initialize map struct
sM = som_map_struct(size(cbv,2), 'msize', [n m], 'rect', 'sheet');
sM.codebook = [];
% convert 
for k=1:m
  sM.codebook = [sM.codebook; cbv(k:m:num_n,:)];
end

% convert hits
n = size(hits,1);
f = zeros(n,m);
for k=1:m
  f(:,k) = hits(:,k);
end
