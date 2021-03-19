% function cl = clhgraph(x,k,sfct)
%
% DESCRIPTION
%   provides cluster labels 1 to k from hypergraph partitioning 
%   sfct is ignored
%
% Copyright (c) 1998-2002 by Alexander Strehl


function cl = clhgraph(x,k,sfct)
% TODO - use 3rd parameter to hmetis - weights
cl = hmetis(x,k);
