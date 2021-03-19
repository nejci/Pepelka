function [W,distances] = compute_relation(data,scale_sig,order)
%
%      [W,distances] = compute_relation(data,scale_sig) 
%       Input: data= Feature_dimension x Num_data
%       ouput: W = pair-wise data similarity matrix
%              Dist = pair-wise Euclidean distance
%
%
% Jianbo Shi, 1997 

% For 2D data
% N = size(data,2);
% distances = zeros(N,N);
% for j = 1:N
%   distances(j,:) = (sqrt((data(1,:)-data(1,j)).^2 +...
%                 (data(2,:)-data(2,j)).^2));
% end

% MUST use sqrt or the results are not so good.
distances = sqrt(X2distances(data'));

if ~exist('scale_sig','var') || isempty(scale_sig) 
    scale_sig = 0.05*max(distances(:));
end

if ~exist('order','var')
  order = 2;
end

tmp = (distances/scale_sig).^order;

W = exp(-tmp);

