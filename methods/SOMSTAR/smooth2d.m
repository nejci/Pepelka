function Ys = smooth2d(Y,theta)
% Smooth values in Y using 2D kernel.
% It is Matlab implementation (simplified) of R function "smooth.2d" from
% package "fields".
% We assume that values in Y come from a grid with a step of 1 in both
% directions.
%
% Description from R docs: 
% "An approximate Nadaraya Watson kernel smoother is obtained by first
% discretizing the locations to a grid and then using convolutions to find
% and to apply the kernel weights. The main advantage of this function is a
% smoother that avoids explicit looping."
%
% Y:     vector or matrix 
% theta: scalar or vector of length size(Y,2); every column of Y will be 
%        divided by theta.
%
% Nejc Ilc, 2014

[xdim,ydim] = size(Y);

if ~exist('theta','var') || isempty(theta)
   theta = 1; 
end
if length(theta) ~= 2
   theta = [theta,theta]; 
end

NN = ones(xdim,ydim);
m = xdim;
n = ydim;
M = 2*m;
N = 2*n;

xg = combvec(1:M,1:N)';
center = [M/2,N/2];

% scale coordinates by theta and compute the kernel
x1 = bsxfun(@rdivide,xg,theta);
x2 = bsxfun(@rdivide,center,theta);
rdist = sqdistance2(x1,x2);
out = reshape(exp(-rdist),M,N);

temp = zeros(M,N);
temp(center(1),center(2)) = 1;
wght = fftn(out)./ (fftn(temp)*M*N);

temp = zeros(M,N);
temp(1:m,1:n) = Y;
temp(isnan(temp)) = 0;

temp = fftn(temp) .* wght;
temp2 = real(ifftn(temp)*M*N);
temp2 = temp2(1:m,1:n);

temp = zeros(M,N);
temp(1:m,1:n) = NN;
temp(isnan(temp)) = 0;
temp = fftn(temp) .* wght;
temp3 = real(ifftn(temp)*M*N);
temp3 = temp3(1:m,1:n);

Ys = temp2./temp3;
