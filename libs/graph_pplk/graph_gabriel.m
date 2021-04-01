% GABRIEL: Finds the Gabriel connectivity graph among points in P dimensions.
%
%     Usage: [G,connect,dist] = graph_gabriel(data,{noplot},{tol})
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
%           connect = [n x n] boolean adjacency matrix.
%           dist =    corresponding edge lengths (Euclidean distances);
%                       non-connected edge distances are given as zero.
%
%
% Gabriel, KR & RR Sokal. 1969. A new statistical approach to geographic variation
%   analysis.  Syst. Zool. 18:259-278.
% Matula, DW & RR Sokal. 1980. Properties of Gabriel graphs relevant to geographic
%   variation research and the clustering of points in the plane.  Geogr. Analysis
%   12:205-222.

% RE Strauss, 1/8/99
%   9/7/99 -   changed plot colors for Matlab v5.
%   10/23/02 - remove restriction of p==2;
%              produce plot only for p==2.
%   6/4/03 -   added tolerance to the minimum A-B circle.

function [G,dist] = graph_gabriel(data,tol)

if (nargin < 2)
    tol = [];
end

if (isempty(tol))
    tol = 1e-10;
end

[n,dim] = size(data);

if (n<2)
    error('GABRIEL: requires at least two points ');
end

dist = sqdistance2(data);

if dim < 4
    % for 2D and 3D cases employ Delaunay triangulation to speed things up.
    % Gabriel is a subgraph of delaunay triangulation.
    DT = DelaunayTri(data);
    E = DT.edges;
    numE = size(E,1);
    d = zeros(n);
    % test all edges whether they belong to Gabriel graph
    for e = 1:numE
        i = E(e,1);
        j = E(e,2);
        dij = dist(i,j);
        d(i,j) = dij;
        % is there a point k that lies in disk i-j?
        for k = 1:n  
            if (k~=i && k~=j)
                % if such k exists, there is no edge between i and j
                if dist(i,k)+dist(j,k)-tol <= dij
                    d(i,j) = 0; 
                    break;
                end
            end
        end
    end
    
else
    d = triu(dist);
    % Cycle thru all possible pairs of points
    for i = 1:(n-1)
        for j = (i+1):n
            dij = dist(i,j);
            for k = 1:n
                if (k~=i && k~=j)
                    if (dist(i,k)+dist(j,k)-tol <= dij)
                        d(i,j) = 0;
                        break;
                    end
                end
            end
        end
    end
end
% form graph format
G = sparse(tril(d'));
