function umat = computeUmat(map,smoothing)

if ~exist('smoothing','var') || isempty(smoothing)
   smoothing = 0; 
end

codebook = map.codebook;
d = squareform(pdist(codebook));
x = map.topol.msize(1);
y = map.topol.msize(2);

umat = zeros(x,y);

if x<2 || y <2
    error('Cannot compute Umat for map with dimension < 2.');
end

% convert (x,y) indices to linear indices
L = reshape(1:x*y,x,y);

% check if the map is larger than 2 x 2 (otherwise it is only corners)
if x > 2 || y > 2
    % inner nodes
    for ix = 2:(x-1)
        for iy = 2:(y-1)
            S = ...
                d(L(ix,iy),L(ix-1,iy-1)) + ...
                d(L(ix,iy),L(ix,iy-1)) + ...
                d(L(ix,iy),L(ix+1,iy-1)) + ...
                d(L(ix,iy),L(ix+1,iy)) + ...
                d(L(ix,iy),L(ix+1,iy+1)) + ...
                d(L(ix,iy),L(ix,iy+1)) + ...
                d(L(ix,iy),L(ix-1,iy+1)) + ...
                d(L(ix,iy),L(ix-1,iy));
            umat(ix,iy) = S/8;
        end
    end
    % left border
    for ix = 2:(x-1)
        iy = 1;
        S = ...
            d(L(ix,iy),L(ix+1,iy)) + ...
            d(L(ix,iy),L(ix+1,iy+1)) + ...
            d(L(ix,iy),L(ix,iy+1)) + ...
            d(L(ix,iy),L(ix-1,iy+1)) + ...
            d(L(ix,iy),L(ix-1,iy));
        umat(ix,iy) = S/5;
    end
    % right border
    for ix = 2:(x-1)
        iy = y;
        S = ...
            d(L(ix,iy),L(ix-1,iy-1)) + ...
            d(L(ix,iy),L(ix,iy-1)) + ...
            d(L(ix,iy),L(ix+1,iy-1)) + ...
            d(L(ix,iy),L(ix+1,iy)) + ...
            d(L(ix,iy),L(ix-1,iy));
        umat(ix,iy) = S/5;
    end
    % top border
    for iy = 2:(y-1)
        ix = 1;
        S = ...
            d(L(ix,iy),L(ix,iy-1)) + ...
            d(L(ix,iy),L(ix+1,iy-1)) + ...
            d(L(ix,iy),L(ix+1,iy)) + ...
            d(L(ix,iy),L(ix+1,iy+1)) + ...
            d(L(ix,iy),L(ix,iy+1));
        umat(ix,iy) = S/5;
    end
    % bottom border
    for iy = 2:(y-1)
        ix = x;
        S = ...
            d(L(ix,iy),L(ix-1,iy-1)) + ...
            d(L(ix,iy),L(ix,iy-1)) + ...
            d(L(ix,iy),L(ix,iy+1)) + ...
            d(L(ix,iy),L(ix-1,iy+1)) + ...
            d(L(ix,iy),L(ix-1,iy));
        umat(ix,iy) = S/5;
    end
end

% compute umat values for corners
if x >= 2 && y >= 2
    % top left corner
    ix = 1;
    iy = 1;
    S = ...
        d(L(ix,iy),L(ix+1,iy)) + ...
        d(L(ix,iy),L(ix+1,iy+1)) + ...
        d(L(ix,iy),L(ix,iy+1));
    umat(ix,iy) = S/3;
    
    % bottom left corner
    ix = x;
    iy = 1;
    S = ...
        d(L(ix,iy),L(ix,iy+1)) + ...
        d(L(ix,iy),L(ix-1,iy+1)) + ...
        d(L(ix,iy),L(ix-1,iy));
    umat(ix,iy) = S/3;
    
    % top right corner
    ix = 1;
    iy = y;
    S = ...
        d(L(ix,iy),L(ix,iy-1)) + ...
        d(L(ix,iy),L(ix+1,iy-1)) + ...
        d(L(ix,iy),L(ix+1,iy));
    umat(ix,iy) = S/3;
    
    % bottom right corner
    ix = x;
    iy = y;
    S = ...
        d(L(ix,iy),L(ix-1,iy-1)) + ...
        d(L(ix,iy),L(ix,iy-1)) + ...
        d(L(ix,iy),L(ix-1,iy));
    umat(ix,iy) = S/3;
end

% smooth the Umat
if smoothing > 0
    umat = smooth2d(umat,smoothing);
elseif smoothing < 0
    error('Bad value for smoothing parameter');
end






