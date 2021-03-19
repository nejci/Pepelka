function DI=dunns2(data,labels,distance)   
%%%Dunn's index for clustering compactness and separation measurement
% dunns2(data,labels,distance)
% data = data matrix [n X d] 
% labels = [n X 1] vector of data labels
% distance  =   string that is passed to pdist to calculate distance between
%               data (e.g. 'euclidean', 'correlation', ...). If [],
%               'euclidean' is used as default.
% -------------------------------------------------------------------------
%
% Copyright (c) 2010, Julian Ramos
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

if ~exist('distance','var') || isempty(distance)
    distance = 'euclidean';
end

i= length(unique(labels));
distM = squareform(pdist(data,distance));
ind = labels;

denominator=[];

for i2=1:i
    indi=find(ind==i2);
    indj=find(ind~=i2);
    x=indi;
    y=indj;
    temp=distM(x,y);
    denominator=[denominator;temp(:)];
end

num=min(min(denominator)); 
neg_obs=zeros(size(distM,1),size(distM,2));

for ix=1:i
    indxs=find(ind==ix);
    neg_obs(indxs,indxs)=1;
end

dem=neg_obs.*distM;
dem=max(max(dem));

DI=num/dem;
end