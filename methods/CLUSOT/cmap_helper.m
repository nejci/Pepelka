function [cmap]=cmap_helper(name,n)
% [cmap]=cmap_helper(name,n)
%
% cmap_helper creates colormap name
% with size n
%
% Input:
% name - colormap name as string
% n    - size of colormap
%
% Output:
% cmap - colormap
% 
% $Id: cmap_helper.m 369 2005-02-04 11:40:18Z dome $
% D. Brugger, 29 January 2005
% util/cmap_helper.m

%if(nargin == 0)
%  test_cmap_helper();
%  return
%end

switch name
 case 'autumn'
  cmap = autumn(n);
 case 'bone'
  cmap = bone(n);
 case 'colorcube'
  cmap = colorcube(n);
 case 'cool'
  cmap = cool(n);
 case 'copper'
  cmap = copper(n);
 case 'flag'
  cmap = flag(n);
 case 'gray'
  cmap = gray(n);
 case 'hot'
  cmap = hot(n);
 case 'hsv'
  cmap = hsv(n);
 case 'jet'
  cmap = jet(n);
 case 'lines'
  cmap = lines(n);
 case 'pink'
  cmap = pink(n);
 case 'prism'
  cmap = prism(n);
 case 'spring'
  cmap = spring(n);
 case 'summer'
  cmap = summer(n);
 case 'white'
  cmap = white(n);
 otherwise
  cmap = winter(n);
end

function [t]=test_cmap_helper()

names={'autumn','bone','colorcube','cool','copper','flag','gray','hot', ...
       'hsv','jet','lines','pink','prism','spring','summer', ...
       'white','winter'};
N=10;
for k=1:length(names)
  cmap = cmap_helper(names{k},N);
  if(sum(size(cmap) == [N,3]) ~= 2)
    error('!!! test_cmap_helper failed !!!')
  end
end

fprintf('test_cmap_helper succeded\n')