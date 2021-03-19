function [Q,R,K,fident] = parseinputs(varargin)
% Check input and output
% Copyright (c) 2009, Yi Cao
narginchk(1,3);

Q=varargin{1};

if nargin<2
    R=Q;
    fident = true;
else
    fident = false;
    R=varargin{2};
end

if isempty(R)
    fident = true;
    R=Q;
end

if ~fident
    fident = isequal(Q,R);
end

if nargin<3
    K=1;
else
    K=varargin{3};
end