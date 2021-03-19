% function cl = cspa(cls,k)
%
% DESCRIPTION
%  Performs CSPA for CLUSTER ENSEMBLES
%
% Copyright (c) 1998-2002 by Alexander Strehl

function cl = cspa(cls,k,weights)

N = size(cls,2);

if ~exist('k','var') || isempty(k)
    k = max(max(cls));
end

weightsModeON = 1;
if ~exist('weights','var') || isempty(weights)
    weightsModeON = 0;
end

if weightsModeON
    s = computeCO(cls',[],'labels','full',weights);
    
    % Normalize similarities on [0,1].
    % Diagonal is 0, thus MIN is 0 and this normalization does not produce
    % additional zeros.
    MAX = max(s(:));
    MIN = min(s(:));    
    if MAX > 1 || MIN < 0
        s = (s-MIN) ./ (MAX-MIN);
    end
    s(1:N+1:N^2) = 1;    
    s = checks(s);
else
    clbs = clstoclbs(cls);
    s = clbs' * clbs;
    s = checks(s./size(cls,1));
end

cl = metis(s,k);

