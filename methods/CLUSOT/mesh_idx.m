function pos=mesh_idx(np,res)
% [pos]=mesh_idx(np,res)
%
% mesh_idx computes the mesh index for neurons at positions np,
% where res is the resolution of the mesh
%
% Inputs:
% np - (matrix) neuron positions (n x 2)-matrix
% res - (double) resolution of mesh
%
% Output:
% pos - (matrix) mesh positions (n x 2)-matrix
% 
% $Id: mesh_idx.m 456 2005-03-20 11:38:58Z dome $
% D. Brugger, 20 March 2005
% util/mesh_idx.m

% compute mesh idx for neurons at positions np
pos=round(np.*(1/res)+1);

