function [A,B] = compute_relation(data,scale_sig,order)
%
%      [W,Dist] = compute_relation(data,scale_sig) 
%       Input: data= Feature_dimension x Num_data
%       ouput: W = pair-wise data similarity matrix
%              Dist = pair-wise Euclidean distance
%
%
% Jianbo Shi, 1997 

if (~exist('order')),
  order = 2;
end

B = zeros(length(data),length(data));
for j = 1:length(data),
  B(j,:) = (sqrt((data(1,:)-data(1,j)).^2 +...
                (data(2,:)-data(2,j)).^2));
end

if (~exist('scale_sig')),
    scale_sig = 0.05*max(B(:));
end

tmp = (B/scale_sig).^order;

A = exp(-tmp);

