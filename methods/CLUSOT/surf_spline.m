function [f,dfdx,dfdy]=surf_spline(np,p,hn,ix,iy,res,bflag)
% [f,dfdx,dfdy]=surf_spline(np,p,hn,ix,iy,res,bflag)
%
% surf_spline computes the surface function and partial
% derivatives of a modified gaussian, where the variance
% is influenced by a cubic spline through control points p.
%
% Inputs:
% np - (vector) neuron position of central neuron, (1 x 2) vector
% p - (matrix) control points, (n x 2) matrix
% hn - (double) normalized frequency of neuron np
% ix - (vector) x dimension of euclidian plane, [x_low,x_high]
% iy - (vector) y dimension of euclidian plane, [y_low,y_high]
% res - (double) resolution used for surface 
% bflag - (bool) if set to true the handling of special cases
%                 uses mirrored neuron positions. Default: true [optional]
%
% Outputs:
% sf - (matrix) surface function
% dfdx, dfdy - (matrix) partial derivatives
%
% $Id: surf_spline.m 578 2005-04-18 09:51:22Z dome $
% D. Brugger, 14 March 2005
% algo/surf_spline.m

if(nargin < 7)
  bflag=true;
end

% prepare grids
[X,Y]=meshgrid(ix(1):res:ix(2),iy(1):res:iy(2));
% compute angles of all points in plane...
angles=atan2(Y-np(2),X-np(1));
% ... and to neighboring units
a=atan2(p(:,2)-np(2),p(:,1)-np(1))';

% bring angles to [0,2*pi]
angles=angles+pi;
a=a+pi;

% sort angles into increasing order
[sa,idx]=sort(a);
% sort p's likewise
pp=p(idx,:);

% check for special cases
for k=1:size(sa,2)
  if(k == size(sa,2))
    alpha = (sa(1)+2*pi)-sa(k);
  else
    alpha = sa(k+1)-sa(k);
  end
  
  % handle special cases
  if(~bflag) % (a) traditional way
    if(alpha > pi)
      % fprintf('Handling special case ''alpha > pi''...');
      alpha_np=mod( alpha/2+sa(k), 2*pi );
      % insert position of np and knot alpha_np
      sa = [sa alpha_np];
      pp = [pp; np];
      % sort again
      [sa,idx]=sort(sa);
      pp=pp(idx,:);
    else
      if(alpha == pi && size(p,1) == 2)
        % fprintf('Handling special case ''alpha == pi && size(p,1) == 2''...');
        % insert np twice
        alpha_np1=sa(k) + (pi/2);
        alpha_np2=mod((sa(k+1) + (pi/2)),2*pi);
        sa = [sa(k) sa(k+1) alpha_np1 alpha_np2];
        pp = [pp; np; np];
        % sort again
        [sa,idx]=sort(sa);
        pp=pp(idx,:);
      end
    end
  else % (b) using mirrored neurons
    % indices of points alpha=sa(i2)-sa(i1)
    if(k == size(sa,2))
      i2=1; i1=k;
    else
      i2=k+1; i1=k;
    end
    if(alpha > pi)
%      fprintf('Handling special case ''alpha > pi'' ...')
      alpha_i1=mod(sa(i1)+pi,2*pi);
      alpha_i2=mod(sa(i2)+pi,2*pi);
      sa = [sa alpha_i1 alpha_i2];
      pp = [pp; (pp(i1,:)-np)*-1+np; (pp(i2,:)-np)*-1+np];
      % sort again
      [sa,idx]=sort(sa);
      pp=pp(idx,:);
    elseif(alpha == pi)
      if(size(p,1) > 2)
%        fprintf(['Handling special case ''alpha == pi && size(p,1) ' ...
%                 '> 2'' ...'])
        % index of next point in order, if needed wrap around
        if(i2 == size(sa,2))
          i3 = 1;
        else
          i3 = i2+1;
        end
        alpha_i3=mod(sa(i3)+pi,2*pi);
        sa = [sa alpha_i3];
        pp = [pp; (pp(i3,:)-np)*-1+np];
        % sort again
        [sa,idx]=sort(sa);
        pp=pp(idx,:);
      elseif(size(p,1) == 2)
%        fprintf(['Handling special case ''alpha == pi && size(p,1) ' ...
%                 ' == 2'' ...'])
        d1=edist(pp(1,:),np);
        d2=edist(pp(2,:),np);
        % compute position of "phantom neuron"
        if(d1 < d2)
          alpha_ph = mod(sa(1)+pi/2,2*pi);
          ph = rotate2d(pp(1,:)-np,pi/2)+np;        
        else
          alpha_ph = mod(sa(2)+pi/2,2*pi);
          ph = rotate2d(pp(2,:)-np,pi/2)+np;
        end
        alpha_ph2 = mod(alpha_ph+pi,2*pi);
        sa = [sa alpha_ph alpha_ph2];
        pp = [pp; ph; (ph-np)*-1+np];
        % sort again
        [sa,idx]=sort(sa);
        pp=pp(idx,:);
      end
    end
  end
end

% add last knot position and point as we compute 
% a closed curve and adjust angles to cover [0,2*pi)
alpha_min=sa(1); 
sa = sa-alpha_min;
sa=[sa 2*pi];
angles=angles-alpha_min;
[r,c]=find(angles < 0);
for k=1:size(r,1)
  angles(r(k),c(k))=angles(r(k),c(k))+2*pi;
end

% make colorplot w/r to angle
%figure;
%surf(angles);
%xyz_label();
%title('Angles');
pp=[pp; pp(1,:)];
[a,b,c,d]=splines(pp,sa);

% plot computed spline
%plot_spline(a,b,c,d,pp,sa,[sa(1):0.01:sa(end)],'Spline',0,np)

% TODO: Modify spline interpolation function to handle whole matrices
% now compute sigma
sigma=zeros(size(angles));
%size(angles)
pts=spline_ipol(sa,a,b,c,d,angles);
%size(pts)
sigma=sqrt((pts(:,:,1)-np(1)).^2 + (pts(:,:,2)-np(2)).^2);
%size(sigma)
%Ax=zeros(size(angles)); Ay=zeros(size(angles));
%Bx=zeros(size(angles)); By=zeros(size(angles));
%Cx=zeros(size(angles)); Cy=zeros(size(angles));
%Dx=zeros(size(angles)); Dy=zeros(size(angles));
%Sx=zeros(size(angles)); Sy=zeros(size(angles));
%TMPx=zeros(size(angles)); TMPy=zeros(size(angles));
%for k=1:size(angles,1)
%  [pts,au,bu,cu,du,tmp1]=spline_ipol(sa,a,b,c,d,angles(k,:));
%  Ax(k,:)=au(:,1)'; Ay(k,:)=au(:,2)';
%  Bx(k,:)=bu(:,1)'; By(k,:)=bu(:,2)';
%  Cx(k,:)=cu(:,1)'; Cy(k,:)=cu(:,2)';
%  Dx(k,:)=du(:,1)'; Dy(k,:)=du(:,2)';
%  TMPx(k,:)=tmp1(:,1)'; TMPy(k,:)=tmp1(:,2)';
%  Sx(k,:)=pts(:,1)'; Sy(k,:)=pts(:,2)'; 
%  sigma(k,:)=sqrt((pts(:,1)-np(1)).^2 + (pts(:,2)-np(2)).^2)';
%end

% handle special case 'alpha == pi && size(p,1) == 2'
[r,c]=find(sigma == 0);
for k=1:size(r,1)
  sigma(r(k),c(k))=eps;
end

% compute gradient
sigma2 = sigma.^2;
%sigma4 = sigma2.^2;
%iX = np(1) - X; jY = np(2) - Y;
Xi = X - np(1); Yj = Y - np(2); 

f=hn.*exp(-(Xi.^2+Yj.^2)./(2*sigma2));
%figure; surf(f); title('Control function'); xyz_label();

% for now! compute gradient numerically

[dfdx,dfdy]=gradient(f);

% TODO: have to deal with some numeric issues, when
% calculating the gradient analytically
%[Fx,Fy]=gradient(atan2(Yj,Xi));

%dfdx=zeros(size(f));
%dfdx = ...
%    f.*((iX ./ sigma2) + ...
%        (((Xi.^2 + Yj.^2) .* -1 .* (np(1) - Sx + np(2) - Sx)  .* ...
%          ((3.*Dx.*(TMPx.^2)) + (2.*Cx.*TMPx) + Bx) .*  Fx) ./ sigma4));
%dfdy = zeros(size(f));
%dfdy = ...
%    f.*((jY ./ sigma2) + ...
%        (((Xi.^2 + Yj.^2) .* -1 .* (np(1) - Sy + np(2) - Sy) .* ...
%          ((3.*Dy.*(TMPy.^2)) + (2.*Cy.*TMPy) + By) .* Fy) ./ sigma4));