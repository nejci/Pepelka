function [B]=cluster_border(M,T,ix,iy,slope,g)
% [B]=cluster_border(M,T,ix,iy,slope,g)
%
% TODO: <short description>
%
% TODO: <explain vars>
% 
% $Id: cluster_border.m 476 2005-03-23 11:34:59Z dome $
% D. Brugger, 28 February 2005
% algo/cluster_border.m

if(nargin == 0)
  test_cluster_border();
  return
end

v=[0:pi/T:pi-(pi/T)];
%U=zeros(length(v),2); L=zeros(length(v),2);
U=[]; L=[]; pos_u=1; pos_l=1;
for k=1:length(v)
  [s1,s2]=linseg(M,v(k),ix,iy);

  if(~isempty(s1))
    r = bresenham(M,s1);
    l_set=0;
    for l=1:size(r,1)
      x=r(l,1); y=r(l,2);
      if(slope(y,x) < g)
        L(k,:)=r(l,:);
        l_set=1;
      else
        % point (x,y) with slope(y,x) >= g
        % does *not* belong to cluster
        break;
      end
    end
    if(l_set)
      pos_l = pos_l + 1;
    end
  end

  if(~isempty(s2))
    r = bresenham(M,s2); found_flag=0;
    u_set=0;
    for l=1:size(r,1)
      x=r(l,1); y=r(l,2);
      if(slope(y,x) < g)
        U(pos_u,:)=r(l,:);
        u_set=1;
      else
        % point (x,y) with slope(y,x) >= g
        % does *not* belong to cluster
        break;
      end
    end
    if(u_set)
      pos_u = pos_u + 1;
    end
  end
end

% b contains first "upper" border points
% followed by "lower" border points
B = [U; L];

% "border" does at least include local maximum
if(isempty(B))
  B = M;
end

function test_cluster_border()
sf=peaks(30); ix=[1,30]; iy=[1,30];
mask=imregionalmax(sf);
M=sf.*mask; [r,c]=find(M > 0)
[gx,gy]=gradient(sf);
slope=sqrt(gx.^2+gy.^2);

% visualize results
figure;
z=-10; cmap=hsv(size(r,1));
g=1.0; T=0;
for ll=1:5
subplot(3,2,ll);
surf(sf); hold on;
T = 2^ll;
for k=1:size(r,1)
  B=cluster_border([r(k), c(k)], T, ix, iy, slope, g);
  scatter3(r(k), c(k), z, 30, cmap(k,:), 'Marker', 'o');
  for l=1:size(B,1)
    scatter3(B(l,1),B(l,2),z, 30, cmap(k,:), 'Marker', 'x');
  end
  plot3([B(:,1); B(1,1)],[B(:,2); B(1,2)],ones(size(B,1)+1,1)*z,'Color',cmap(k,:))
end
title(sprintf('Cluster borders for T=%d',T));
xlabel('x'); ylabel('y'); zlabel('z');
end

subplot(3,2,6);
surfc(slope);
title('Slope');
xlabel('x'); ylabel('y'); zlabel('z');
