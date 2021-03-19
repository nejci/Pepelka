function [conn,conn_norm]=indexCON(data,labels,L)
% [conn,conn_norm,conn_inv]=indexCON(data,labels,L)
%--------------------------------------------------------------------------
% Cluster internal validity index - Connectivity index
% Index should be minimized.
%--------------------------------------------------------------------------
% INPUTS
%   data		(matrix)	matrix [n X d] with n d-dimensional samples
%
%	labels		(vector)	array of non-negative integers determining
%							the labels of data samples
%	L			(scalar)    is a parameter determining the number of
%                           neighbours of every data point that contribute
%                           to the connectivity measure
%--------------------------------------------------------------------------
% OUTPUTS:
%   conn         (scalar)	value of the Connectivity index
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2013,  Nejc Ilc
%
%------- REFERENCE --------------------------------------------------------
% Handl, J., & Knowles, J. (2005). Exploiting the Trade-off --The Benefits
% of Multiple Objectives in Data Clustering. In C. Coello Coello, A.
% Hernandez Aguirre, & E. Zitzler (Eds.), Evolutionary Multi-Criterion
% Optimization (Vol. 3410, pp. 547�560). Springer Berlin / Heidelberg.
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 22-July-2013 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <nejc.ilc@fri.uni-lj.si>
%==========================================================================


[N,D]=size(data);

% 1. construct L-nearest neighbour graph on data
Lnn = knnsearch(data,[],L);
% Lnn is a matrix [nPatterns x L], where each i-th row contains the list of
% indices of the L-closest data points to i

% conn=0;
% for i=1:N
%     %cluster label of i
%     l_i=labels(i);
%     
%     for j=1:L
%         l_j=labels(Lnn(i,j));
%         if l_i ~= l_j
%             conn=conn + 1/j;
%         end
%     end
% end

% Fully vectorized version (FAST)
LAB=labels(Lnn);
diff = bsxfun(@ne,labels,LAB);
conn = sum(sum(bsxfun(@rdivide,diff,1:L),2));


maxConn=N*sum(1./(1:L));

% Normalized conn index on interval [0,1]
conn_norm=conn/maxConn;

% %inverse conn index, as used by Vega-Pons 2010
% if conn==0
%     conn_inv=1; % to avoid division by zero
% else
%     conn_inv=1/conn;
% end