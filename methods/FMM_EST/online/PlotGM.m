function h=PlotGM(mixEM,meansEM,varsEM,varargin)
% function h=PlotGM(mix,means,vars,threshold,nPoints,R)
%
% Input:
%  Mixture:
%   mix                  - 1 x nModes matrix
%   means                - nDimensions x nModes matrix
%   vars                 - nDimensions x nDimensions x nModes array
%  Parameters:
%   threshold (optional) - plot only significant components mix(iMode)>threshold (default: threshold=0)
%   nPoints (optional)   - number of points in the graph (default: nPoints=100)
%   R (optional)         - nDimensions x nDimensionsNew matrix for mapping to a lower-dimensional
%                          nDimensionsNew space (R'*data). If not supplied and nDimensions>3 first
%                          2 components are taken.
%TO DO:
%   type(optional)       - 1(default) - r standard deviations line (2D and 3D only)
%                          2 - full pdf (1D and 2D only)
%                          3 - r standard deviations largest eigen vector (2D and 3D only)
%                          4 - r standard deviations all eigen vectors (2D and 3D only)
%   r(optional)          - plot previous r standard deviations from the means (default 2)
%
% Output:
%  h                    - graphics handle
%
% Author:   Z.Zivkovic
% Date:     10-2-2003

nD=size(meansEM,1);
nModes=length(mixEM);

%defaults
nPoints=100;
fThreshold=0;
if (nargin>3)
    fThreshold=varargin{1};
elseif (nargin>4)
    nPoints=varargin{2};
end

%reduce dimensionality if needed
if (nD>3)
    if (nargin>5)
        R=varargin{3};
        if (~isempty(R))
            %map
            nD=size(R,2);
            meansEM=R'*meansEM;
            varsEMnew=zeros(nD,nD,nModes);
            for iModes=1:nModes
                varsEMnew(:,:,iModes)=R'*varsEM(:,:,iModes)*R;
            end
            varsEM=varsEMnew;
        end
    end
end
if (nD>3)%check again
    %map
    nD=2;
    meansEM=meansEM(1:nD,:);
    varsEM=varsEM(1:nD,1:nD,:);
end

%plot
switch nD
    case 1
        %pdf
        stdEM=sqrt(squeeze(varsEM));
        
        %borders
        if (size(stdEM,1)~=size(meansEM,1))
            stdEM=stdEM';
        end
        fStart=min(meansEM-3*stdEM);
        fStop=max(meansEM+3*stdEM);
        fStep=(fStop-fStart)/(nPoints-1);
        dataEM=fStart:fStep:fStop;
        
        %generate points
        rpdfPi=zeros(nModes,nPoints);
        rpdfP=zeros(1,nPoints);
        for iData=1:nPoints
            pdfPi=zeros(1,nModes);
            pdfPtot=0;
            diff2=(meansEM-dataEM(iData)).^2;
            for iModes=1:nModes
                if (mixEM(iModes)>fThreshold)
                    %gauss
                    dist=0.5*diff2(iModes)/varsEM(iModes);
                    pdfPi(iModes)=mixEM(iModes)*(1/(sqrt(2*pi)*stdEM(iModes)))*exp(-dist);
                    pdfPtot=pdfPtot+pdfPi(iModes);
                end
            end
            rpdfPi(:,iData)=pdfPi';
            rpdfP(iData)=pdfPtot;
        end
        h=plot(dataEM,rpdfP,'k-');
        %set(h,'LineWidth',2,'Color',[0 0 0]);
        %plot modes
        hold on
        for iModes=1:nModes
            if (mixEM(iModes)>fThreshold)
                plot(dataEM,rpdfPi(iModes,:),'k:')
            end
        end
        %hold off
        
    case 2
        %2-sigma contours
        w=0:2*pi/nPoints:2*pi;
        rCircle=2*[cos(w);sin(w)];%circle
        nPoints=length(w);
        hold on;
        for iModes=1:nModes
            if (mixEM(iModes)>fThreshold)
                [Fi,Lambda]=eig(varsEM(:,:,iModes));
                K=Fi*Lambda^0.5;
                rPlot=K*rCircle+meansEM(:,iModes)*ones(1,nPoints);
                plot(rPlot(1,:),rPlot(2,:),'k-');
                %set(h,'LineWidth',3)
            end
        end
        hold off;
        
    case 3
        %plot main axis
        hold on;
        for iModes=1:nModes
            if (mixEM(iModes)>fThreshold)
                [Fi,Lambda]=eig(varsEM(:,:,iModes));
                if(0)
                    %not sorted !!!! find max!!!
                    riMax=find(max(diag(Lambda))==diag(Lambda));
                    iMax=riMax(1);
                    vecMax=1*Fi(:,iMax)*Lambda(iMax,iMax)^0.5;
                    Line3D=[meansEM(:,iModes)-vecMax meansEM(:,iModes)+vecMax];
                    h=plot3(Line3D(1,:),Line3D(2,:),Line3D(3,:),'k-');
                    set(h,'LineWidth',3);%,'Color',[0 0 0]);
                else
                    K=Fi*(2*Lambda^0.5);
                    [x,y,z]=sphere(10);
                    nP=size(x,1);
                    rData=K*[reshape(x,1,nP^2);reshape(y,1,nP^2);reshape(z,1,nP^2)]+meansEM(:,iModes)*ones(1,nP^2);
                    h=surf(reshape(rData(1,:),nP,nP),reshape(rData(2,:),nP,nP),reshape(rData(3,:),nP,nP));
                    set(h,'FaceAlpha',0.2,'LineStyle','none','FaceColor',[1 0 0],'FaceLighting','phong');
                end
            end
        end
        hold off;
        
end



