function [idx,D]=knnsearch(varargin)
% KNNSEARCH   Linear k-nearest neighbor (KNN) search
% IDX = knnsearch(Q,R,K) searches the reference data set R (n x d array
% representing n points in a d-dimensional space) to find the k-nearest
% neighbors of each query point represented by each row of Q (m x d array).
% The results are stored in the (m x K) index array, IDX.
%
% IDX = knnsearch(Q,R) takes the default value K=1.
%
% IDX = knnsearch(Q) or IDX = knnsearch(Q,[],K) does the search for R = Q.
%
% Rationality
% Linear KNN search is the simplest appraoch of KNN. The search is based on
% calculation of all distances. Therefore, it is normally believed only
% suitable for small data sets. However, other advanced approaches, such as
% kd-tree and delaunary become inefficient when d is large comparing to the
% number of data points. On the other hand, the linear search in MATLAB is
% relatively insensitive to d due to the vectorization. In  this code, the
% efficiency of linear search is further improved by using the JIT
% aceeleration of MATLAB. Numerical example shows that its performance is
% comparable with kd-tree algorithm in mex.
%
% See also, kdtree, nnsearch, delaunary, dsearch

% By Yi Cao at Cranfield University on 25 March 2008
%--------------------------------------------------------------------------
%Copyright (c) 2009, Yi Cao
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%--------------------------------------------------------------------------

% Check inputs
[Q,R,K,fident] = parseinputs(varargin{:});

% Check outputs
error(nargoutchk(0,2,nargout));

% C2 = sum(C.*C,2)';
[N,M] = size(Q);
L=size(R,1);
idx = zeros(N,K);
D = idx;

if K==1
    % Loop for each query point
    for k=1:N
        d=zeros(L,1);
        for t=1:M
            d=d+(R(:,t)-Q(k,t)).^2;
        end
        if fident
            d(k)=inf;
        end
        [D(k),idx(k)]=min(d);
    end
else
    for k=1:N
        d=zeros(L,1);
        for t=1:M
            d=d+(R(:,t)-Q(k,t)).^2;
        end
        if fident
            d(k)=inf;
        end
        [s,t]=sort(d);
        idx(k,:)=t(1:K);
        D(k,:)=s(1:K);
    end
end
if nargout>1
    D=sqrt(D);
end
