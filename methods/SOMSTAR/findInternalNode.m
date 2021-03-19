function [bestx,besty,xcoords,ycoords] = findInternalNode(umat,ix,iy,xcoords,ycoords)

[xdim,ydim] = size(umat);

% first we check if the current position is already associated
% with a centroid.  if so, simply return the coordinates of that centroid
if xcoords(ix,iy)>-1 && ycoords(ix,iy)>-1
    bestx = xcoords(ix,iy);
    besty = ycoords(ix,iy);
    return;
end

% try to find a smaller value in the immediate neighborhood make our
% current position the square with the minimum value. if a minimum value
% other that our own current value cannot be found then we are at a
% minimum. search the neighborhood; three different cases: inner element,
% corner element, side element
minVal = umat(ix,iy);
minX = ix;
minY = iy;

% (ix,iy) is an inner map element
if (ix > 1 && ix < xdim && iy > 1 && iy < ydim)
    if (umat(ix-1,iy-1) < minVal)
        minVal = umat(ix-1,iy-1);
        minX = ix-1;
        minY = iy-1;
    end
    if (umat(ix,iy-1) < minVal)
        minVal = umat(ix,iy-1);
        minX = ix;
        minY = iy-1;
    end
    if (umat(ix+1,iy-1) < minVal)
        minVal = umat(ix+1,iy-1);
        minX = ix+1;
        minY = iy-1;
    end
    if (umat(ix+1,iy) < minVal)
        minVal = umat(ix+1,iy);
        minX = ix+1;
        minY = iy;
    end
    if (umat(ix+1,iy+1) < minVal)
        minVal = umat(ix+1,iy+1);
        minX = ix+1;
        minY = iy+1;
    end
    if (umat(ix,iy+1) < minVal)
        minVal = umat(ix,iy+1);
        minX = ix;
        minY = iy+1;
    end
    if (umat(ix-1,iy+1) < minVal)
        minVal = umat(ix-1,iy+1);
        minX = ix-1;
        minY = iy+1;
    end
    if (umat(ix-1,iy) < minVal)
        minVal = umat(ix-1,iy);
        minX = ix-1;
        minY = iy;
    end
    
    % (ix,iy) is top left corner
elseif (ix == 1 && iy == 1)
    if (umat(ix+1,iy) < minVal)
        minVal = umat(ix+1,iy);
        minX = ix+1;
        minY = iy;
    end
    if (umat(ix+1,iy+1) < minVal)
        minVal = umat(ix+1,iy+1);
        minX = ix+1;
        minY = iy+1;
    end
    if (umat(ix,iy+1) < minVal)
        minVal = umat(ix,iy+1);
        minX = ix;
        minY = iy+1;
    end
    
    % (ix,iy) is bottom left corner
elseif (ix == xdim && iy == 1)
    if (umat(ix,iy+1) < minVal)
        minVal = umat(ix,iy+1);
        minX = ix;
        minY = iy+1;
    end
    if (umat(ix-1,iy+1) < minVal)
        minVal = umat(ix-1,iy+1);
        minX = ix-1;
        minY = iy+1;
    end
    if (umat(ix-1,iy) < minVal)
        minVal = umat(ix-1,iy);
        minX = ix-1;
        minY = iy;
    end
    
    % (ix,iy) is bottom right corner
elseif (ix == xdim && iy == ydim)
    if (umat(ix-1,iy-1) < minVal)
        minVal = umat(ix-1,iy-1);
        minX = ix-1;
        minY = iy-1;
    end
    if (umat(ix,iy-1) < minVal)
        minVal = umat(ix,iy-1);
        minX = ix;
        minY = iy-1;
    end
    if (umat(ix-1,iy) < minVal)
        minVal = umat(ix-1,iy);
        minX = ix-1;
        minY = iy;
    end
    
    % (ix,iy) is top right corner
elseif (ix == 1 && iy == ydim)
    if (umat(ix,iy-1) < minVal)
        minVal = umat(ix,iy-1);
        minX = ix;
        minY = iy-1;
    end
    if (umat(ix+1,iy-1) < minVal)
        minVal = umat(ix+1,iy-1);
        minX = ix+1;
        minY = iy-1;
    end
    if (umat(ix+1,iy) < minVal)
        minVal = umat(ix+1,iy);
        minX = ix+1;
        minY = iy;
    end
    
    % (ix,iy) is a top side element
elseif (ix == 1  && iy > 1 && iy < ydim)
    if (umat(ix,iy-1) < minVal)
        minVal = umat(ix,iy-1);
        minX = ix;
        minY = iy-1;
    end
    if (umat(ix+1,iy-1) < minVal)
        minVal = umat(ix+1,iy-1);
        minX = ix+1;
        minY = iy-1;
    end
    if (umat(ix+1,iy) < minVal)
        minVal = umat(ix+1,iy);
        minX = ix+1;
        minY = iy;
    end
    if (umat(ix+1,iy+1) < minVal)
        minVal = umat(ix+1,iy+1);
        minX = ix+1;
        minY = iy+1;
    end
    if (umat(ix,iy+1) < minVal)
        minVal = umat(ix,iy+1);
        minX = ix;
        minY = iy+1;
    end
    
    % (ix,iy) is a left side element
elseif (ix > 1 && ix < xdim && iy == 1 )
    if (umat(ix+1,iy) < minVal)
        minVal = umat(ix+1,iy);
        minX = ix+1;
        minY = iy;
    end
    if (umat(ix+1,iy+1) < minVal)
        minVal = umat(ix+1,iy+1);
        minX = ix+1;
        minY = iy+1;
    end
    if (umat(ix,iy+1) < minVal)
        minVal = umat(ix,iy+1);
        minX = ix;
        minY = iy+1;
    end
    if (umat(ix-1,iy+1) < minVal)
        minVal = umat(ix-1,iy+1);
        minX = ix-1;
        minY = iy+1;
    end
    if (umat(ix-1,iy) < minVal)
        minVal = umat(ix-1,iy);
        minX = ix-1;
        minY = iy;
    end
    
    % (ix,iy) is a bottom side element
elseif (ix == xdim && iy > 1 && iy < ydim)
    if (umat(ix-1,iy-1) < minVal)
        minVal = umat(ix-1,iy-1);
        minX = ix-1;
        minY = iy-1;
    end
    if (umat(ix,iy-1) < minVal)
        minVal = umat(ix,iy-1);
        minX = ix;
        minY = iy-1;
    end
    if (umat(ix,iy+1) < minVal)
        minVal = umat(ix,iy+1);
        minX = ix;
        minY = iy+1;
    end
    if (umat(ix-1,iy+1) < minVal)
        minVal = umat(ix-1,iy+1);
        minX = ix-1;
        minY = iy+1;
    end
    if (umat(ix-1,iy) < minVal)
        minVal = umat(ix-1,iy);
        minX = ix-1;
        minY = iy;
    end
    
    
    % (ix,iy) is a right side element
elseif (ix > 1 && ix < xdim && iy == ydim)
    if (umat(ix-1,iy-1) < minVal)
        minVal = umat(ix-1,iy-1);
        minX = ix-1;
        minY = iy-1;
    end
    if (umat(ix,iy-1) < minVal)
        minVal = umat(ix,iy-1);
        minX = ix;
        minY = iy-1;
    end
    if (umat(ix+1,iy-1) < minVal)
        minVal = umat(ix+1,iy-1);
        minX = ix+1;
        minY = iy-1;
    end
    if (umat(ix+1,iy) < minVal)
        minVal = umat(ix+1,iy);
        minX = ix+1;
        minY = iy;
    end
    if (umat(ix-1,iy) < minVal)
        minVal = umat(ix-1,iy);
        minX = ix-1;
        minY = iy;
    end
end

%if successful
% move to the square with the smaller value, i.e., call find.internal.node on this new square
% note the RETURNED x-y coords in the x.coords and y.coords matrix at the current location
% return the RETURNED x-y coordinates
if (minX ~= ix || minY ~= iy)    
    [bestx,besty,xcoords,ycoords] = findInternalNode(umat,minX,minY,xcoords,ycoords);
    %if (explicit)
    % if explicit is set show the exact connected component
    % otherwise construct a connected componenent where all
    % nodes are connected to a central node
    %	x.coords[ix,iy] <<- min.x
    %	y.coords[ix,iy] <<- min.y
    %	list(bestx=min.x,besty=min.y)
    %
    %else
    xcoords(ix,iy) = bestx;
    ycoords(ix,iy) = besty;
    return;
    %end
   
else
    % we have found a minimum
    % note the current x-y in the xcoords and ycoords matrix
    % return the current x-y coordinates
    xcoords(ix,iy) = ix;
    ycoords(ix,iy) = iy;
    bestx=ix;
    besty=iy;
    return;
end

