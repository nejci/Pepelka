function [G,dist] = graph_gabriel(data,tol,dist)
% GRAPH_GABRIEL: Finds the Gabriel connectivity graph among N points in P dimensions.
%
%     Usage: [G,connect,dist] = graph_gabriel(data,{tol},{dist})
%
%           data =    [N x P] matrix of point coordinates.
%           tol  =    optional tolerance for squared distance from a third point
%                     to the minimum A-B circle (i.e., that circle on whose
%                     circumference A & B are at opposite points) [default = 1e-6].
%           dist =    [N x N] matrix of precomputed pairwise point distances
%           -------------------------------------------------------------------------
%			G =		  [n x n] sparse matrix that represents a graph. Nonzero
%					  entries in matrix G represent the weights of the edges.
%           dist =    corresponding edge lengths (Euclidean distances);
%                     non-connected edge distances are given as zero.
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
%   6/4/03   - added tolerance to the minimum A-B circle.
% Nejc Ilc
%   11/19/13 - faster implementation
%            - added parameter dist - a precomuputed distance matrix
%   9/18/15  - added handling of data with 1 or 2 points

if (nargin < 2)
    tol = [];
    dist = [];
end

if (isempty(tol))
    tol = 1e-10;
end

[n,dim] = size(data);

if isempty(dist)
    dist = sqdistance2(data);
end


if n == 1
    % only one data point
    %error('GABRIEL: requires at least two points ');
    G = sparse(1,1); % zero scalar    
    
elseif n == 2
    % two points
    G = sparse(2,1,dist(2,1),2,2); % distance between two points in lower left corner
    
else
    % More than two points
    if dim < 4
        % for 2D and 3D cases employ Delaunay triangulation to speed things up.
        % Gabriel is a subgraph of Delaunay triangulation.
        DT = delaunayTriangulation(data);
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
        % MEX implementation
        d = gabriel(dist,tol);
        
        % Pure Matlab code - slow
        % d = triu(dist);
        % % Cycle thru all possible pairs of points
        % for i = 1:(n-1)
        %     for j = (i+1):n
        %         dij = dist(i,j);
        %         for k = 1:n
        %             if (k~=i && k~=j)
        %                 if (dist(i,k)+dist(j,k)-tol <= dij)
        %                     d(i,j) = 0;
        %                     break;
        %                 end
        %             end
        %         end
        %     end
        % end
    end
    
    % form a graph format (diagonal elements are 0)
    G = sparse(tril(d',-1));
end
