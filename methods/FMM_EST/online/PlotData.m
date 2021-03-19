function h=PlotData(dataEM,varargin)
% function h=PlotData(data,R)
%
% Input:
%
%  data          - nDimensions x nSamples matrix
%  R (optional)  - nDimensions x nDimensionsNew matrix for mapping to a lower-dimensional
%                  nDimensionsNew space (R'*data). If not supplied and nDimensions>3 first 
%                  2 components are taken.
%  type(optional)- standard matlab plot type - see help plot (default 'b.')
%
% Output:
%  h             - graphics handle 
%
% Author:   Z.Z.
% Date:     10-2-2003

nD=size(dataEM,1);

%reduce dimensionality if needed
if (nD>3)
    if (nargin>1)
        R=varargin{1};
        if (~isempty(R))
            %map
            dataEM=R'*dataEM;
            nD=size(R,2);
        end
    end
end
    
if (nD>3)%check again
  %map
  nD=2;        
  dataEM=dataEM(1:nD,:);
end

%plot type
type='.';
if (nargin==3)
    type=varargin{2};
end

%plot
switch nD
case 1
    %histogram
    %borders
    nData=size(dataEM,2);
    nPoints=50;%parameter !!!
    fStart=min(dataEM);
    fStop=max(dataEM);
    fStep=(fStop-fStart)/(nPoints-1);
    edgesEM=fStart:fStep:fStop;
    HistEM=histc(dataEM,edgesEM);
    h=bar(edgesEM,HistEM/(fStep*nData));
    set(h,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[1 1 1]);
   
case 2
    h=plot(dataEM(1,:),dataEM(2,:),type);
    set(h,'Color',[0.6 0.6 0.6]);
    
case 3
    h=plot3(dataEM(1,:),dataEM(2,:),dataEM(3,:),type);
    set(h,'Color',[0.6 0.6 0.6],'MarkerSize',5);
    grid on
    
end

        
    
    