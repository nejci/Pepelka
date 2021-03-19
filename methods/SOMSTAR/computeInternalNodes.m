function coords = computeInternalNodes(umat)

[xdim,ydim] = size(umat);
xcoords = ones(xdim,ydim)*-1;
ycoords = ones(xdim,ydim)*-1;
maxVal = max(umat(:));

% iterate over the map and find the centroid for each element
for ix = 1:xdim
    for iy = 1:ydim
        [~,~,xcoords,ycoords] = findInternalNode(umat,ix,iy,xcoords,ycoords);
    end
end

coords.xcoords = xcoords;
coords.ycoords = ycoords;