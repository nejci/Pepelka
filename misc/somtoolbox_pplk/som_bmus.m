function [Bmus,Qerrors] = som_bmus(sMap, sData, which_bmus, mask)
% Works for correlation distance also! Distance function should be
% specified in sMap.trainhist(2).distfun
%
%SOM_BMUS Find the best-matching units from the map for the given vectors.
%
% [Bmus, Qerrors] = som_bmus(sMap, sData, [which], [mask])
% 
%   bmus = som_bmus(sM,sD);
%   [bmus,qerrs] = som_bmus(sM,D,[1 2 3]);
%   bmus = som_bmus(sM,D,1,[1 1 0 0 1]);
%
%  Input and output arguments ([]'s are optional): 
%   sMap     (struct) map struct
%            (matrix) codebook matrix, size munits x dim
%   sData    (struct) data struct
%            (matrix) data matrix, size dlen x dim
%   [which]  (vector) which BMUs are returned, [1] by default 
%            (string) 'all', 'best' or 'worst' meaning [1:munits],
%                     [1] and [munits] respectively  
%   [mask]   (vector) mask vector, length=dim, sMap.mask by default
%
%   Bmus     (matrix) the requested BMUs for each data vector, 
%                     size dlen x length(which)
%   Qerrors  (matrix) the corresponding quantization errors, size as Bmus
%
% NOTE: for a vector with all components NaN's, bmu=NaN and qerror=NaN
% NOTE: the mask also effects the quantization errors
%
% For more help, try 'type som_bmus' or check out online documentation.
% See also  SOM_QUALITY.

%{
%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% som_bmus
%
% PURPOSE
%
% Finds Best-Matching Units (BMUs) for given data vector from a given map.
%
% SYNTAX
%
%  Bmus = som_bmus(sMap, sData)
%  Bmus = som_bmus(..., which)
%  Bmus = som_bmus(..., which, mask)
%  [Bmus, Qerrs] = som_bmus(...)
%
% DESCRIPTION
%
% Returns the indexes and corresponding quantization errors of the
% vectors in sMap that best matched the vectors in sData.
%
% By default only the index of the best matching unit (/vector) is
% returned, but the 'which' argument can be used to get others as
% well. For example it might be desirable to get also second- and
% third-best matching units as well (which = [1:3]). 
%
% A mask can be used to weight the search process. The mask is used to
% weight the influence of components in the distance calculation, as
% follows: 
%
%   distance(x,y) = (x-y)' diag(mask) (x-y)
%
% where x and y are two vectors, and diag(mask) is a diagonal matrix with 
% the elements of mask vector on the diagonal. 
%
% The vectors in the data set (sData) can contain unknown components
% (NaNs), but the map (sMap) cannot. If there are completely empty
% vectors (all NaNs), the returned BMUs and quantization errors for those 
% vectors are NaNs.
%
% REQUIRED INPUT ARGUMENTS
%
%   sMap              The vectors from among which the BMUs are searched
%                     for. These must not have any unknown components (NaNs).
%            (struct) map struct
%            (matrix) codebook matrix, size munits x dim
%                     
%   sData             The data vector(s) for which the BMUs are searched.
%            (struct) data struct
%            (matrix) data matrix, size dlen x dim
%
% OPTIONAL INPUT ARGUMENTS 
%
%   which    (vector) which BMUs are returned, 
%                     by default only the best (ie. which = [1])
%            (string) 'all', 'best' or 'worst' meaning [1:munits],
%                     [1] and [munits] respectively  
%   mask     (vector) mask vector to be used in BMU search, 
%                     by default sMap.mask, or ones(dim,1) in case
%                     a matrix was given
%
% OUTPUT ARGUMENTS
% 
%   Bmus     (matrix) the requested BMUs for each data vector, 
%                     size dlen x length(which)
%   Qerrors  (matrix) the corresponding quantization errors, 
%                     size equal to that of Bmus
%
% EXAMPLES
%
% Simplest case:
%  bmu = som_bmus(sM, [0.3 -0.4 1.0]);
%           % 3-dimensional data, returns BMU for vector [0.3 -0.4 1]
%  bmu = som_bmus(sM, [0.3 -0.4 1.0], [3 5]);
%           % as above, except returns the 3rd and 5th BMUs
%  bmu = som_bmus(sM, [0.3 -0.4 1.0], [], [1 0 1]);
%           % as above, except ignores second component in searching
%  [bmus qerrs] = som_bmus(sM, D);
%           % returns BMUs and corresponding quantization errors 
%           % for each vector in D
%  bmus = som_bmus(sM, sD);
%           % returns BMUs for each vector in sD using the mask in sM
%
% SEE ALSO
% 
%  som_quality      Measure the quantization and topographic error of a SOM.
%
% Copyright (c) 1997-2000 by the SOM toolbox programming team.
% http://www.cis.hut.fi/projects/somtoolbox/
%
% Version 1.0beta juuso 071197, 101297 
% Version 2.0alpha juuso 201198 080200
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%}
% check arguments and initialize

error(nargchk(1, 4, nargin));  % check no. of input args is correct

% sMap
if isstruct(sMap), 
  switch sMap.type, 
   case 'som_map', M = sMap.codebook; 
   case 'som_data', M = sMap.data;
   otherwise, error('Invalid 1st argument.');
  end
else 
  M = sMap; 
end
[munits dim] = size(M);
if any(any(isnan(M))), 
  error ('Map codebook must not have missing components.');
end

isNormalized = 0;
% data
if isstruct(sData), 
  switch sData.type, 
   case 'som_map', D = sData.codebook;
   case 'som_data'
       D = sData.data; 
       if isfield(sData,'isNormalized')
           isNormalized = sData.isNormalized;
       end           
   otherwise, error('Invalid 2nd argument.');
  end
else 
  D = sData;
end
[dlen ddim] = size(D);
if dim ~= ddim, 
  error('Data and map dimensions do not match.')
end

% which_bmus
if nargin < 3 | isempty(which_bmus) | any(isnan(which_bmus)), 
  which_bmus = 1; 
else
  if ischar(which_bmus), 
    switch which_bmus,
     case 'best', which_bmus = 1; 
     case 'worst', which_bmus = munits; 
     case 'all', which_bmus = [1:munits];
    end
  end
end

% mask
if nargin < 4 | isempty(mask) | any(isnan(mask)), 
  if isstruct(sMap) & strcmp(sMap.type,'som_map'), 
    mask = sMap.mask; 
  elseif isstruct(sData) & strcmp(sData.type,'som_map'), 
    mask = sData.mask; 
  else
    mask = ones(dim,1); 
  end
end
if size(mask,1)==1, mask = mask'; end
if all(mask == 0), 
  error('All components masked off. BMU search cannot be done.');
end

% Added by Nejc
% distance - passed as sMap.trainhist(2).distfun
distfun = 'sqEuclidean';
som_mode='euc';
if isstruct(sMap)
	if isfield(sMap,'trainhist')
		if isfield(sMap.trainhist,'distfun')
			if length(sMap.trainhist) >1
				distfun = sMap.trainhist(2).distfun;
			end
		end
		if isfield(sMap.trainhist,'distmode')
			if length(sMap.trainhist) >1
				som_mode = strtok(sMap.trainhist(2).distmode,'_');
			end
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% action

Bmus = zeros(dlen,length(which_bmus));
Qerrors = Bmus;

% The BMU search involves calculating weighted Euclidian distances or
% correlation
% to all map units for each data vector. Basically this is done as
%   for i=1:dlen, 
%     for j=1:munits, 
%       for k=1:dim,
%         Dist(j,i) = Dist(j,i) + mask(k) * (D(i,k) - M(j,k))^2;
%       end
%     end
%   end
% where mask is the weighting vector for distance calculation. However, taking 
% into account that distance between vectors m and v can be expressed as
%   |m - v|^2 = sum_i ((m_i - v_i)^2) = sum_i (m_i^2 + v_i^2 - 2*m_i*v_i)
% this can be made much faster by transforming it to a matrix operation:
%   Dist = (M.^2)*mask*ones(1,d) + ones(m,1)*mask'*(D'.^2) - 2*M*diag(mask)*D'
% 
% In the case where there are unknown components in the data, each data
% vector will have an individual mask vector so that for that unit, the 
% unknown components are not taken into account in distance calculation.
% In addition all NaN's are changed to zeros so that they don't screw up 
% the matrix multiplications.

% calculate distances & bmus

% This is done a block of data at a time rather than in a
% single sweep to save memory consumption. The 'Dist' matrix has 
% size munits*blen which would be HUGE if you did it in a single-sweep
% operation. If you _want_ to use the single-sweep version, just 
% set blen = dlen. If you're having problems with memory, try to 
% set the value of blen lower. 
blen = min(munits,dlen);

% check available memory
[userview systemview] = memory;
% how much do we need (W1 + WD (diag) + dconst + bmus + ddist + Dist)
requiredMemBytes = ddim*dlen*8 + ddim*ddim*8 + dlen*8 + dlen*8 + dlen*8 + munits*dlen*8;
freeBytes = min(systemview.PhysicalMemory.Available, userview.MemAvailableAllArrays);

memSaving = 0;
if freeBytes < requiredMemBytes
	memSaving = 1;
end
memSaving = 1;

% handle unknown components
Known = ~isnan(D);
D(find(~Known)) = 0; 
unknown = find(sum(Known')==0); % completely unknown vectors 

switch distfun
	case {'sqEuclidean', 'default'}
		if ~memSaving
			% constant matrices
			W1 = (mask*ones(1,dlen)) .* Known';	
			WD = 2*diag(mask)*D';   % constant matrix
			dconst = ((D.^2)*mask); % constant term in the distances
		else
			dconst = 0;
		end
		
		i0 = 0; 
		while i0+1<=dlen, 
		  % calculate distances 
		  inds = [(i0+1):min(dlen,i0+blen)]; 
		  i0 = i0+blen;      
		  
		  if memSaving
			  Dist = zeros(munits,length(inds));
			  % do it in a loop to save memory (slower, we do not compute WD beforehand)
			  for i=1:length(inds)
				  for j=1:munits
					  Dist(j,i) = ((D(inds(i),:) - M(j,:)).^2) * mask;
				  end
			  end
			  
		  else
			  Dist = (M.^2)*W1(:,inds) - M*WD(:,inds); % plus dconst for each sample
		  end
		  
		  % find the bmus and the corresponding quantization errors
		  if all(which_bmus==1)
			  [Q B] = min(Dist); 
		  else
			  [Q B] = sort(Dist); 
		  end
		  if munits==1
			  Bmus(inds,:) = 1; 
		  else
			  Bmus(inds,:) = B(which_bmus,:)'; 
		  end
		  
		  if memSaving
			  Qerrors(inds,:) = Q(which_bmus,:)';
		  else
			  Qerrors(inds,:) = Q(which_bmus,:)' + dconst(inds,ones(length(which_bmus),1));
		  end 
		end
		
		
	case 'correlation'
		
		% transforming codebook vectors
		% in a case of dist_mode=='dotprod', M is already normalized
		% so we can skip this phase optimization.
		if ~strcmp(som_mode,'dotprod')
			Mt = M - repmat(mean(M,2),1,dim);
			Mnorm = sqrt(sum(Mt.^2, 2));
			Mt = Mt ./ Mnorm(:,ones(1,dim));
		else 
			Mt = M;
		end
		
		i0 = 0;
		while i0+1<=dlen			
			inds = [(i0+1):min(dlen,i0+blen)];
			i0 = i0+blen;
			% calculate distances
			% transforming portion of data vectors, only if not normalized
			% before (gSOM)
			if ~isNormalized
				Dt = D(inds,:) - repmat(mean(D(inds,:),2),1,ddim);
				Dnorm = sqrt(sum(Dt.^2, 2));
				Dt = Dt ./ Dnorm(:,ones(1,ddim));
			else 
				Dt = D(inds,:);
			end
			
			% calculate Pearson coefficient and transform it to distance
			Dist = (-1)* (Mt * Dt');
			
			% find the bmus and the corresponding quantization errors
			if all(which_bmus==1)
				[Q B] = min(Dist);
			else
				[Q B] = sort(Dist);
			end
			
			if munits==1
				Bmus(inds,:) = 1;
			else
				Bmus(inds,:) = B(which_bmus,:)';
			end
			
			Qerrors(inds,:) = Q(which_bmus,:)'+1; % "correlation" distance on [0,2]
		end
				
		
	case 'ucorrelation'
		% transforming codebook vectors
		% in a case of dist_mode=='dotprod', M is already normalized
		% so we can skip this phase - optimization.
		Mt = M;
		if ~strcmp(som_mode,'dotprod')
			Mnorm = sqrt(sum(Mt.^2, 2));
			Mt = Mt ./ Mnorm(:,ones(1,dim));		
		end
		
		i0 = 0; 
		while i0+1<=dlen, 
		  % calculate distances 
		  inds = [(i0+1):min(dlen,i0+blen)]; 
		  i0 = i0+blen;      
		  
		  % transforming portion of data vectors
		  Dt = D(inds,:);
		  if ~isNormalized
			Dnorm = sqrt(sum(Dt.^2, 2));
			Dt = Dt ./ Dnorm(:,ones(1,ddim));
		  end
		  
		  % calculate dot product
		  Dist = (Mt * Dt');
		  
		  % find the bmus and the corresponding quantization errors
		  if all(which_bmus==1) 
			  [Q B] = max(Dist); 
		  else
			  [Q B] = sort(Dist,'descend'); 
		  end
		  
		  if munits==1 
			  Bmus(inds,:) = 1; 
		  else
			  Bmus(inds,:) = B(which_bmus,:)'; 
		  end
		  
		  Qerrors(inds,:) = (-1)*Q(which_bmus,:)'+1; % "correlation" distance on [0,2]
		end
		
	otherwise 
		error('Wrong distance function specified!');
	
end

 

% completely unknown vectors
if ~isempty(unknown), 
  Bmus(unknown,:) = NaN;
  Qerrors(unknown,:) = NaN;
end

Qerrors = sqrt(Qerrors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
