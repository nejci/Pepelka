function [mixEM,meansEM,varsEM]=EMRandomInit(dataEM,nModes)
% [mixEM,meansEM,varsEM]=EMRandomInit(data,nModes) 
%
% Input:
% - data    - nDimensions x nData vector - some samples
% - nModes  - initial number of modes
%
% Output:
% - random initialization:
%   mix     - 1 x nModes matrix
%   means   - nDimensions x nModes matrix
%   vars    - nDimensions x nDimensions x nModes array
%
% Author:   Z.Z.
% Date:     10-2-2003
nD=size(dataEM,1);
nData=size(dataEM,2);

mixEM=(1/nModes)*ones(1,nModes);meansEM=zeros(nD,nModes);varsEM=zeros(nD,nD,nModes);
var0=(0.1/nD)*trace(cov(dataEM'))*eye(nD);
for iModes=1:nModes
    iData=floor(nData*rand)+1;
    meansEM(:,iModes)=dataEM(:,iData);
    varsEM(:,:,iModes)=var0;
end