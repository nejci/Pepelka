function [C,BB,Mx,pMx,BP]=clusot_cluster(sf,slope,g,T,bflag)
% [C,BB,Mx,pMx,BP]=clusot_cluster(sf,slope,g,T,bflag)
%
% clusot_cluster computes clusters using surface function sf
% and its slope. A cluster border is detected when slope exceeds
% given parameter g. Resolution of the border is influenced by T.
%
% Inputs:
% sf - (matrix) surface function
% slope - (matrix) slope of sf
% g - (double) threshhold parameter
% T - (double) resolution parameter for scanlines
% bflag - (bool) if set to true, clusters are automerged if
%                 their cluster regions overlap in C. Default: false.
% Outputs:
% C - (matrix) cluster matrix
% BB - (matrix) cluster border points (n x 2) matrix
% Mx - (vector) all local maxima
% pMx - (vector) processed local maxima, e.g. those which weren't
%                ignored 
% BP - (matrix) all cluster border points as determined by
%               interpolation method between points in BB (n x 2)
%               matrix
% 
% $Id: clusot_cluster.m 848 2005-07-25 20:39:25Z dome $
% D. Brugger, 01 March 2005
% algo/clusot_cluster.m

if(nargin == 0)
  test_clusot_cluster();
  return;
end

if(nargin < 5)
  bflag=false;
end

[n,m]=size(sf);
ix=[1 m]; iy=[1 n];
mask=imregionalmax(sf);
[r,c]=find(sf.*mask > 0);
% Mx contains x- and y-coordinates
% in the first two columns and the value
% of the local maximum in the third column
Mx=zeros(length(r),3);
for k=1:length(r)
  Mx(k,:)=[c(k) r(k) sf(r(k),c(k))];
end

% sort Mx into descending order (therfore 'flipdim') 
% using value of maximum, e.g. column 3, as key
Mx=flipdim(sortrows(Mx,3),1);

% init cluster matrix and cluster number
C=zeros(n,m); cn=1;
% compute cluster matrix
pos=1; BB={}; BP={};
for k=1:length(Mx)
  % check if local maximum is not already in cluster
  % - if so do nothing
  x = Mx(k,1); y = Mx(k,2);
  if(C(y,x) == 0)
    fprintf('Processing local maximum at (%d,%d)...\n',x,y);
    % compute cluster border
    B = cluster_border(Mx(k,1:2), T, ix, iy, slope, g);
    % remember processed maxima and border
    pMx(cn,:)=Mx(k,:); BB{cn}=B; 
    % interpolate (linear for now)
    B = lin_ipol(B);
    BP{cn}=B; 
    % fill
    if(bflag) % automerge clusters if necessary
      cnum=unique(diag(C(B(:,2),B(:,1))));
      cnum=cnum(find(cnum > 0));
      if(isempty(cnum)) % nothing to do
        C=fill_with_number(C,B,0,cn);
        cn = cn + 1;
      else % automerge
        fprintf('Automergeing clusters...\n');
        [cnum_min,min_idx]=min(cnum);
        C=fill_with_number(C,B,0,cn);
        C(find(C==cn))=cnum_min;
        for l=1:size(cnum,1)
          C(find(C==cnum(l)))=cnum_min;
        end
        cnum_all=[cnum; cn];
        
        % clean up C, this assures that cluster numbers form a
        % sequence, e.g. 1,2,3,4,... instead of 1,3,5,10
        cnum(min_idx)=0; [cnum_max,max_idx]=max(cnum);
        while(cnum_max > 0)
          idxs=find(C>cnum_max);
          C(idxs)=C(idxs)-1;
          cnum(max_idx)=0;
          [cnum_max,max_idx]=max(cnum);
        end

        % clean up BP, this assures that merged clusters have a
        % common boundary
        Bn=fix_boundaries({BP{cnum_all}},C,0,cnum_min);
        % set new border points
        BB_all=BB{cnum_min};
        for l=1:size(cnum_all,1)
          BP{cnum_all(l)}=[]; 
          BB_all = union(BB_all, BB{cnum_all(l)}, 'rows');
          BB{cnum_all(l)}=[];
        end
        BP{cnum_min}=Bn;
        BB{cnum_min}=BB_all;
        pMx=setdiff(pMx, pMx(cn,:),'rows');
        pMx=setdiff(pMx, pMx(setdiff(cnum,max(cnum)),:),'rows');
%        BP
        BP=remove_empty(BP);
%        BB
        BB=remove_empty(BB);
        
        cn = max(max(C)) + 1;
      end
    else
      C=fill_with_number(C,B,0,cn);
      cn = cn + 1;
    end
  else
    fprintf('Ignoring local maximum at (%d,%d)...\n',x,y);
  end
end

function C=fill_with_number(C,B,null,cn)
C=vector_edge_fill(C,B,0,cn);
for k=1:size(B,1)
  C(B(k,2),B(k,1))=cn;
end

% use to fix boundaries after merging of clusters
function Bn=fix_boundaries(Bs,C,null,cn)
P=cell(1,size(Bs,2));
for k=1:size(Bs,2)
  P{k}=all_candidates(Bs{k},C,null,cn);
end
idx=1;
while(~all_empty(P,idx))
  if(~isempty(P{idx}))
    tmp=P{idx};
    Bn{idx}=find_contour(C,tmp(1,:),cn);
%    P_before=P
    P{idx}=[];
    P=remove_candidates(P,idx,Bn{idx});
%    P_after=P
    idx=idx+1;
  else
    idx=idx+1;
  end
end

% remove candidates for find_contour which already belong to
% computed contour K
function Ps=remove_candidates(Ps,idx,K)
for k=idx+1:size(Ps,2)
%  fprintf('Removing candidates...\n')
%  Ps{k}
%  K
  Ps{k} = setdiff(Ps{k},K,'rows');
end

% checks if all sets of points in Ps are empty
function b=all_empty(Ps,idx)
for k=idx:size(Ps,2)
%  Ps{k}
  if(~isempty(Ps{k}))
    b=false; %    fprintf('all_empty = false\n');
    return;
  end
end
b=true;
%fprintf('all_empty = true\n');

% determine all points on border Bi which are candidates for being
% starting points in the find_contour function
function P=all_candidates(Bi,C,null,cn)
pos=1; [n,m]=size(C);
% might be stored in a cell if cluster already participated in a merge
if(iscell(Bi)) 
  for l=1:size(Bi,2)
    Btmp=Bi{l};
    for k=1:size(Btmp,1)
      four_neigh=[Btmp(k,1)-1 Btmp(k,2)];
      if(four_neigh(1) < 1 || ...
         (four_neigh(1) >= 1 && four_neigh(1) <= m && ...
          four_neigh(2) >= 1 && four_neigh(2) <= n && ...
          C(four_neigh(2), four_neigh(1)) == null) ...
        )
        P(pos,:)=Btmp(k,:); pos=pos+1;
      end
    end
  end
else % Bi is a matrix
  for k=1:size(Bi,1)
    four_neigh=[Bi(k,1)-1 Bi(k,2)];
    if(four_neigh(1) < 1 || ...
       (four_neigh(1) >= 1 && four_neigh(1) <= m && ...
        four_neigh(2) >= 1 && four_neigh(2) <= n && ...
        C(four_neigh(2), four_neigh(1)) == null) ...
       )
      P(pos,:)=Bi(k,:); pos=pos+1;
    end
  end
end

function test_clusot_cluster()
sf=peaks(30); 
[gx,gy]=gradient(sf);
slope=sqrt(gx.^2+gy.^2);
figure; surfc(sf); xyz_label();
% plot slope
figure;
surfc(slope);
xlabel('x'); ylabel('y'); zlabel('z'); title('Slope');
% set parameters
g=1.0; T=128; 
[C,BB,Mx,pMx]=clusot_cluster(sf,slope,g,T);
plot_helper(C,BB,Mx,pMx,sf,g,0);
g=1.5; 
[C,BB,Mx,pMx]=clusot_cluster(sf,slope,g,T);
plot_helper(C,BB,Mx,pMx,sf,g,0);
[C,BB,Mx,pMx]=clusot_cluster(sf,slope,g,T,'automerge');
plot_helper(C,BB,Mx,pMx,sf,g,1);

function plot_helper(C,BB,Mx,pMx,sf,g,bflag)
% plot surface + computed clusters
figure; 
surf(sf); 
hold on; 
z=-10;
cmap=hsv(max(max(C)));
[n,m]=size(C);
for k=1:n
  for l=1:m
    if(C(k,l) ~= 0)
      plot3(l,k,z,'MarkerEdgeColor',cmap(C(k,l),:), ...
            'MarkerFaceColor',cmap(C(k,l),:), 'Marker','o')
    end
  end
end
% Mark local maxima for control purposes
for k=1:length(Mx)
  plot3(Mx(k,1),Mx(k,2),z,'MarkerEdgeColor','b','Marker','x','MarkerSize',20)
end
if(bflag)
  automerge='true';
else
  automerge='false';
end
xyz_label(); title(sprintf('Surface + Clusters g=%g, automerge=%s', g, automerge));
figure;
hold on;
% visualize scan lines 
[n,m]=size(C);
for k=1:n
  for l=1:m
    if(C(k,l) ~= 0)
      plot3(l,k,0,'MarkerEdgeColor',cmap(C(k,l),:), ...
           'MarkerFaceColor',cmap(C(k,l),:), 'Marker','o')
    end
  end
end
n=size(pMx,1)
for k=1:n
  B=BB{k};
  for l=1:length(B)
    plot3([pMx(k,1) B(l,1)],[pMx(k,2), B(l,2)],[0 0],'Color','b');
  end
end
xyz_label(); title(sprintf('Scanlines g=%g, automerge=%s', g, automerge));
