function [h]=plot_clusot(sf,sl,C,BP,pMx,varargin)
% [h]=plot_clusot(sf,sl,C,BP,pMx,varargin)
%
% plot_clusot draws the clusot surface and slope with grayscale
% colors and superimposes the cluster structure.
%
% Inputs:
% sf - (matrix) surface function as computed by clusot_surf
% sl - (matrix) slope of sf
% C - (matrix) cluster matrix as computed by clusot_cluster
% BP - (matrix) border points -||-
% pMx - (matrix) processed local maxima  -||-
%
% [argID,   (string) See below. 
%    value] type of value is dependent on argID
% 
% argID can be:
% 'flag' - (string) draw borders or fill cluster regions?
%         Can be set to 'fill' or 'border'. Default: 'border'.
% 'width' - (int) line width used for drawing cluster
%           borders. Default: 2
%
% Output:
% h - handles to plots
%
% $Id: plot_clusot.m 867 2005-07-28 12:09:20Z dome $
% D. Brugger, 14 March 2005
% viz/plot_clusot.m

% default values
width=2;
flag='border';

% parse options
i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i}, 
      % argument IDs
     case 'flag',      i=i+1; flag = varargin{i}; 
     case 'width',    i=i+1; width = varargin{i}; 
     otherwise argok=0; 
    end
  else
    argok=0;
  end
  if ~argok, 
    disp(['(plot_clusot) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end


h=zeros(2,1);
cmap=hsv(max(max(C)));
% plot surface + computed clusters
figure; 
hold on; 
plot_local_maxima(sf,C,pMx,cmap);
if(strcmp(flag,'border'))
  plot_border(sf,BP,pMx,cmap,width);
else
  plot_filled(sf,C,pMx,cmap);
end
colormap(gray);
h_tmp=surf(sf); 
xyz_label();
title('Surface + Clusters');
axis tight;
set(h_tmp, 'EdgeColor', 'none');
h(1)=gcf;

% plot slope + computed clusters
figure;
hold on; 
plot_local_maxima(sl,C,pMx,cmap);
if(strcmp(flag,'border'))
  plot_border(sl,BP,pMx,cmap,width);
else
  plot_filled(sl,C,pMx,cmap);
end
colormap(gray);
h_tmp=surf(sl); 
xyz_label();
title('Slope + Clusters');
axis tight;
set(h_tmp, 'EdgeColor', 'none');
h(2)=gcf;

% Marks processed local maxima in current plot
function plot_local_maxima(sf,C,pMx,cmap)
for k=1:size(pMx,1)
  cn=C(pMx(k,2),pMx(k,1));
  color=cmap(cn,:);
  plot3(pMx(k,1),pMx(k,2),sf(pMx(k,2),pMx(k,1)), ...
        'MarkerEdgeColor', color,'Marker','x','MarkerSize',20);
  text(pMx(k,1),pMx(k,2),sf(pMx(k,2),pMx(k,1)),sprintf('Cluster #%d',cn));
end

function plot_border(sf,BP,pMx,cmap,width)
n=size(pMx,1);
for k=1:n
  B=BP{k};
  if(iscell(B)) % cluster has more than one border
%    fprintf('Cluster has more than one border...\n');
    for l=1:size(B,2)
      plot_border_helper(sf,B{l},cmap(k,:),width);
    end
  else
    plot_border_helper(sf,B,cmap(k,:),width);
  end
end

function plot_border_helper(sf,B,color,width)
tmp_x=[B(1:end-1,1) B(2:end,1); ...
       B(end,1) B(1,1)];
tmp_y=[B(1:end-1,2) B(2:end,2); ...
       B(end,2) B(1,2)];
tmp_z=[diag(sf(B(1:end-1,2),B(1:end-1,1))) ...
       diag(sf(B(2:end,2),B(2:end,1))); ...
       diag(sf(B(end,2),B(end,1))) ...
       diag(sf(B(1,2),B(1,1)))];
plot3(tmp_x,tmp_y,tmp_z,'Color',color,'LineWidth',width);
%for k=1:size(B,1)-1
%  plot3([B(k,1) B(k+1,1)],[B(k,2) B(k+1,2)], ...
%        [sf(B(k,2),B(k,1)) sf(B(k+1,2),B(k+1,1))], 'Color', color);
%end
%k=size(B,1);
%plot3([B(k,1) B(1,1)],[B(k,2) B(1,2)], ...
%      [sf(B(k,2),B(k,1)) sf(B(1,2),B(1,1))], 'Color', color);

function plot_filled(sf,C,pMx,cmap)
n=size(pMx,1);
for k=1:n
  ind=find(C == k);
  [x,y]=ind2sub(size(C),ind);
  z=diag(sf(x,y));
  plot3([y y],[x x],[z z], 'Color', cmap(k,:));
  %'MarkerEdgeColor',cmap(k,:), ...
%        'MarkerFaceColor',cmap(k,:), 'Marker','o');
end
%[n,m]=size(C);
%for k=1:n
%  for l=1:m
%    if(C(k,l) ~= 0)
%      plot3(l,k,sf(k,l),'MarkerEdgeColor',cmap(C(k,l),:), ...
%           'MarkerFaceColor',cmap(C(k,l),:), 'Marker','o')
%    end
%  end
%end
