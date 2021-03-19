function [h]=plot_clust(map,clust,cms)
% [h]=plot_clust(map,clust,cms)
%
% TODO: <short description>
%
% TODO: <explain vars>
% 
% $Id: plot_clust.m 545 2005-04-06 12:13:11Z dome $
% D. Brugger, 20 February 2005
% viz/plot_clust.m

figure;
ah=axes;
% calculate neuron positions in grid
tmp = som_vis_coords(map.topol.lattice, map.topol.msize);
msize = map.topol.msize;
m=msize(2);n=msize(1);
if(strcmp(map.topol.lattice,'rect'))
  pos = 1;
  for k=1:m
    positions(pos:pos+n-1,:)=flipdim(tmp(pos:pos+n-1,:),1);
    pos = pos + n;
  end
else
  error('!!! Unkown topology ''%s''. Need to fix positions!!!', map.topol.lattice); 
end

[n1,m1]=size(positions);
set(ah, 'XLim', [min(positions(:,1))-1 max(positions(:,1))+1], ...
        'YLim', [min(positions(:,2))-1 max(positions(:,2))+1], ...
        'Visible', 'off');
%        'XTick', [], ...
%        'YTick', []);
for k=1:n1
  text(positions(k,1),positions(k,2),num2str(clust(k)));
end


num_clust=max(clust);
if(nargin < 4)
  cmap = colormap(jet(num_clust));
else
  cmap = cmap_helper(cms,num_clust);
  colormap(cmap);
end
spacer=0.05;
delta_x = 0.5-spacer; delta_y = 0.5-spacer; pos = 1;
for k=1:m
  for l=1:n
    lu = positions(pos,:) + [-delta_x delta_y];
    ll = positions(pos,:) + [-delta_x -delta_y];
    ru = positions(pos,:) + [delta_x delta_y];
    rl = positions(pos,:) + [delta_x -delta_y];
    has_rn=0; has_ln=0; has_un=0; has_dn=0;
    if(rneigh_other_clust(pos,l,n,clust))
      has_rn=1;
    end
    if(lneigh_other_clust(pos,l,n,clust))
      has_ln=1;
    end
    if(uneigh_other_clust(pos,l,n,clust))
      has_un=1;
    end
    if(dneigh_other_clust(pos,l,n,clust))
      has_dn=1;
    end
    
    if(has_rn)
      ru_x=ru(1); rl_x=rl(1);
      ru_y=ru(2); rl_y=rl(2);
      if(~has_un)
        ru_y = ru_y + 2*spacer;
      end
      if(~has_dn);
        rl_y = rl_y - 2*spacer;
      end
      lh=line([ru_x rl_x],[ru_y rl_y]);
      set(lh, 'Color', cmap(clust(pos),:), 'LineWidth', 2);
    end
    
    if(has_ln)
      lu_x=lu(1); ll_x=ll(1);
      lu_y=lu(2); ll_y=ll(2);
      if(~has_un)
        lu_y = lu_y + 2*spacer;
      end
      if(~has_dn)
        ll_y = ll_y - 2*spacer;
      end
      lh=line([lu_x ll_x],[lu_y ll_y]);
      set(lh, 'Color', cmap(clust(pos),:), 'LineWidth', 2);
    end

    if(has_un)
      ru_x=ru(1); lu_x=lu(1);
      ru_y=ru(2); lu_y=lu(2);
      if(~has_rn)
        ru_x=ru_x+2*spacer;
      end
      if(~has_ln)
        lu_x=lu_x-2*spacer;
      end
      lh=line([lu_x ru_x],[lu_y ru_y]);
      set(lh, 'Color', cmap(clust(pos),:), 'LineWidth', 2);
    end

    if(has_dn)
      ll_x=ll(1); rl_x=rl(1);
      ll_y=ll(2); rl_y=rl(2);
      if(~has_rn)
        rl_x = rl_x+2*spacer;
      end
      if(~has_ln)
        ll_x = ll_x-2*spacer;
      end
      lh=line([ll_x rl_x],[ll_y rl_y]);
      set(lh, 'Color', cmap(clust(pos),:), 'LineWidth', 2);
    end
    pos = pos + 1;
  end
end

h=gcf;

% following functions determine if neighbors are in a different cluster
function bool=rneigh_other_clust(pos,l,n,clust)
bool=0; 
% k == m <=> pos+n > length(clust)
if (pos+n <= length(clust) && clust(pos) ~= clust(pos+n)) || ...
      pos+n > length(clust)
  bool=1;
end
  
function bool=lneigh_other_clust(pos,l,n,clust)
bool=0;
% k == 1 <=> pos-n <= 0
if (pos-n > 0 && clust(pos) ~= clust(pos-n)) || ...
      pos-n <= 0
  bool=1;
end

function bool=uneigh_other_clust(pos,l,n,clust)
bool=0;
if (l ~= 1 && clust(pos) ~= clust(pos-1)) || ...
         l == 1
  bool=1;
end

function bool=dneigh_other_clust(pos,l,n,clust)
bool=0;
if (l ~= n && clust(pos) ~= clust(pos+1)) || ...
         l == n
  bool=1;
end

         
