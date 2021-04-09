% RNG: Finds the Relative Neighborhood Graph among points in P dimensions.
%
%     Usage: [G,connect,dist] = graph_rng(data,{noplot},{tol})
%
%           data =    [n x p] matrix of point coordinates.
%           noplot =  optional boolean flag indicating that plot should be
%                       suppressed [default = 0].  Plot of the Gabriel graph
%                       is produced only for p=2.
%           tol =     optional tolerance for squared distance from a third point
%                       to the minimum A-B circle (i.e., that circle on whose
%                       circumference A & B are at opposite points) [default = 1e-6].
%           -------------------------------------------------------------------------
%			G =		  [n x n] sparse matrix that represents a graph. Nonzero
%					  entries in matrix G represent the weights of the edges.
%           dist =    corresponding edge lengths (Euclidean distances);
%                     non-connected edge distances are given as zero.
%

% RE Strauss, 1/8/99
%   9/7/99 -   changed plot colors for Matlab v5.
%   10/23/02 - remove restriction of p==2;
%              produce plot only for p==2.
%   6/4/03 -   added tolerance to the minimum A-B circle.
% N Ilc
%	gabriel -> rng
%   speed up by Delaunay triangulation for 2D and 3D cases.
%   9/18/15 - special cases when number of points is 1 or 2

function [G,dist] = graph_rng(data,tol)
if (nargin < 2)
    tol = [];
end

if (isempty(tol))
    tol = 1e-6;
end


[n,dim] = size(data);

dist = sqdistance2(data);

if n == 1
    % one point
    G = sparse(1,1); % zero scalar
    
elseif n == 2
    % two points
    G = sparse(2,1,dist(2,1),2,2); % distance between two points in lower left corner
    
else
    
    % PAIRS = nchoose2(1:n);
    % IND = sub2ind(size(d),PAIRS(:,1),PAIRS(:,2));
    % DIJ = d(IND);
    
    if dim < 4
        % for 2D and 3D cases employ Delaunay triangulation to speed things up.
        % RNG is a subgraph of delaunay triangulation.
        DT = delaunayTriangulation(data);
        E = DT.edges;
        numE = size(E,1);
        d = zeros(n);
        % check all edges whether they belong to RNG graph
        for e = 1:numE
            i = E(e,1);
            j = E(e,2);
            dij = dist(i,j);
            d(i,j) = dij;
            % is there a point k that is closer to the i (or j) than j (or i)?
            for k = 1:n
                if (k~=i && k~=j)
                    % if such k exists, there is no edge between i and j
                    if max(dist(i,k),dist(j,k))-tol <= dij
                        d(i,j) = 0;
                        break;
                    end
                end
            end
        end
        
    else
        % For higher dimensions employ general algorithm (slow for big N).
        % Cycle thru all possible pairs of points
        d = triu(dist);
        for i = 1:(n-1)
            for j = (i+1):n
                dij = dist(i,j);
                % is there a point k that is closer to the i (or j) than j (or i)?
                for k = 1:n
                    if (k~=i && k~=j)
                        if max(dist(i,k),dist(j,k))-tol <= dij
                            % if such k exists, there is no edge between i and j
                            d(i,j) = 0;
                            break;
                        end
                    end
                end
            end
        end
    end
    
    % form graph format
    G=sparse(tril(d'));
end

