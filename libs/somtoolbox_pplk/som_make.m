function sMap = som_make(D, varargin)

%SOM_MAKE Create, initialize and train Self-Organizing Map.
%
% sMap = som_make(D, [[argID,] value, ...])
%
%  sMap = som_make(D);
%  sMap = som_make(D, 'munits', 20);
%  sMap = som_make(D, 'munits', 20, 'hexa', 'sheet');
%  sMap = som_make(D, 'msize', [4 6 7], 'lattice', 'rect');
%
%  Input and output arguments ([]'s are optional): 
%   D        (matrix) training data, size dlen x dim
%            (struct) data struct
%   [argID,  (string) See below. The values which are unambiguous can 
%    value]  (varies) be given without the preceeding argID.
%
%   sMap     (struct) map struct
%
% Here are the valid argument IDs and corresponding values. The values 
% which are unambiguous (marked with '*') can be given without the
% preceeding argID.
%   'init'       *(string) initialization: 'randinit' or 'lininit' (default)
%   'algorithm'  *(string) training: 'seq' or 'batch' (default) or 'sompak'
%   'munits'      (scalar) the preferred number of map units
%   'msize'       (vector) map grid size
%   'mapsize'    *(string) do you want a 'small', 'normal' or 'big' map
%                          Any explicit settings of munits or msize override this.
%   'lattice'    *(string) map lattice, 'hexa' or 'rect'
%   'shape'      *(string) map shape, 'sheet', 'cyl' or 'toroid'
%   'neigh'      *(string) neighborhood function, 'gaussian', 'cutgauss',
%                          'ep' or 'bubble'
%   'topol'      *(struct) topology struct
%   'som_topol','sTopol' = 'topol'
%   'mask'        (vector) BMU search mask, size dim x 1
%   'name'        (string) map name
%   'comp_names'  (string array / cellstr) component names, size dim x 1
%   'tracking'    (scalar) how much to report, default = 1
%   'training'    (string) 'short', 'default', 'long'
%                 (vector) size 1 x 2, first length of rough training in epochs, 
%                          and then length of finetuning in epochs
%
% For more help, try 'type som_make' or check out online documentation.
% See also SOM_MAP_STRUCT, SOM_TOPOL_STRUCT, SOM_TRAIN_STRUCT,
%          SOM_RANDINIT, SOM_LININIT, SOM_SEQTRAIN, SOM_BATCHTRAIN.          

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% som_make
%
% PURPOSE
%
% Creates, initializes and trains a SOM using default parameters.
%
% SYNTAX
%
%  sMap = som_make(D);
%  sMap = som_make(...,'argID',value,...);
%  sMap = som_make(...,value,...);
%
% DESCRIPTION
%
% Creates, initializes and trains a SOM with default parameters. Uses functions
% SOM_TOPOL_STRUCT, SOM_TRAIN_STRUCT, SOM_DATA_STRUCT and SOM_MAP_STRUCT to come
% up with the default values.
%
% First, the number of map units is determined. Unless they are
% explicitly defined, function SOM_TOPOL_STRUCT is used to determine this.
% It uses a heuristic formula of 'munits = 5*dlen^0.54321'. The 'mapsize'
% argument influences the final number of map units: a 'big' map has 
% x4 the default number of map units and a 'small' map has x0.25 the
% default number of map units. 
%
% After the number of map units has been determined, the map size is 
% determined. Basically, the two biggest eigenvalues of the training
% data are calculated and the ratio between sidelengths of the map grid
% is set to this ratio. The actual sidelengths are then set so that 
% their product is as close to the desired number of map units as
% possible.
%
% Then the SOM is initialized. First, linear initialization along two
% greatest eigenvectors is tried, but if this can't be done (the
% eigenvectors cannot be calculated), random initialization is used
% instead.  After initialization, the SOM is trained in two phases:
% first rough training and then fine-tuning. If the 'tracking'
% argument is greater than zero, the average quantization error and
% topographic error of the final map are calculated.
%
% REQUIRED INPUT ARGUMENTS
%
%  D           The data to use in the training.
%     (struct) A data struct. If a struct is given, '.comp_names' field as 
%              well as '.comp_norm' field is copied to the map struct.
%     (matrix) A data matrix, size dlen x dim. The data matrix may
%              contain unknown values, indicated by NaNs. 
%  
% OPTIONAL INPUT ARGUMENTS 
%
%  argID (string) Argument identifier string (see below).
%  value (varies) Value for the argument (see below).
%
% Here are the valid argument IDs and corresponding values. The values 
% which are unambiguous (marked with '*') can be given without the
% preceeding argID.
%   'init'       *(string) initialization: 'randinit' or 'lininit' (default)
%   'algorithm'  *(string) training: 'seq' or 'batch' (default) or 'sompak'
%   'munits'      (scalar) the preferred number of map units
%   'msize'       (vector) map grid size
%   'mapsize'    *(string) do you want a 'small', 'normal' or 'big' map
%                          Any explicit settings of munits or msize override this.
%   'lattice'    *(string) map lattice, 'hexa' or 'rect'
%   'shape'      *(string) map shape, 'sheet', 'cyl' or 'toroid'
%   'neigh'      *(string) neighborhood function, 'gaussian', 'cutgauss',
%                          'ep' or 'bubble'
%   'topol'      *(struct) topology struct
%   'som_topol','sTopol' = 'topol'
%   'mask'        (vector) BMU search mask, size dim x 1
%   'name'        (string) map name
%   'comp_names'  (string array / cellstr) component names, size dim x 1
%   'tracking'    (scalar) how much to report, default = 1
%   'training'    (string) 'short', 'default' or 'long'
%                 (vector) size 1 x 2, first length of rough training in epochs, 
%                          and then length of finetuning in epochs
%
% OUTPUT ARGUMENTS
% 
%  sMap (struct) the trained map struct
%
% EXAMPLES
%
%  To simply train a map with default parameters: 
%
%   sMap = som_make(D); 
%  
%  With the optional arguments, the initialization and training can be
%  influenced. To change map size, use 'msize', 'munits' or 'mapsize'
%  arguments:  
%
%   sMap = som_make(D,'mapsize','big'); or sMap=som_make(D,'big');
%   sMap = som_make(D,'munits', 100);
%   sMap = som_make(D,'msize', [20 10]); 
%
%  Argument 'algorithm' can be used to switch between 'seq' and 'batch'
%  algorithms. 'batch' is the default, so to use 'seq' algorithm: 
%
%   sMap = som_make(D,'algorithm','seq'); or sMap = som_make(D,'seq'); 
%
%  The 'tracking' argument can be used to control the amout of reporting
%  during training. The argument is used in this function, and it is
%  passed to the training functions. To make the function work silently
%  set it to 0.
%
%   sMap = som_make(D,'tracking',0); 
%
% SEE ALSO
% 
%  som_map_struct   Create a map struct.
%  som_topol_struct Default values for SOM topology.
%  som_train_struct Default values for SOM training parameters.
%  som_randinint    Random initialization algorithm.
%  som_lininit      Linear initialization algorithm.
%  som_seqtrain     Sequential training algorithm.
%  som_batchtrain   Batch training algorithm.

% Copyright (c) 1999-2000 by the SOM toolbox programming team.
% http://www.cis.hut.fi/projects/somtoolbox/
% Version 2.1 nejci (personal branch)
% Version 2.0beta juuso 111199

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check arguments

% D
if isstruct(D) 
  data_name = D.name; 
  comp_names = D.comp_names;
  comp_norm = D.comp_norm;
  %dodal Nejc Ilc
  pca_vec=D.pca_vec;
  pca_val=D.pca_val;
  
  D = D.data;
else 
  data_name = inputname(1);
  sDummy = som_data_struct(D(1,:)); 
  comp_names = sDummy.comp_names;
  comp_norm = sDummy.comp_norm;
end
[dlen dim] = size(D);

% defaults
mapsize = '';
sM = som_map_struct(dim); 
sTopol = sM.topol;
munits = prod(sTopol.msize); % should be zero
mask = sM.mask; 
name = sM.name; 
neigh = sM.neigh; 
tracking = 1;
algorithm = 'batch';
distfun = 'default'; %sqEuclidean
dist_mode = 'euc_pure';
som_mode = 'euc';%euc
initalg = 'lininit';
training = 'default';
sample_order_type='random';
advanced = []; % for training {'name',value}

% varargin
i=1; 
while i<=length(varargin), 
  argok = 1; 
  if ischar(varargin{i}), 
    switch varargin{i}, 
      % argument IDs
     case 'mask',       i=i+1; mask = varargin{i}; 
     case 'munits',     i=i+1; munits = varargin{i}; 
     case 'msize',      i=i+1; sTopol.msize = varargin{i}; 
                        munits = prod(sTopol.msize); 
     case 'mapsize',    i=i+1; mapsize = varargin{i}; 
     case 'name',       i=i+1; name = varargin{i};
     case 'comp_names', i=i+1; comp_names = varargin{i}; 
     case 'lattice',    i=i+1; sTopol.lattice = varargin{i};
     case 'shape',      i=i+1; sTopol.shape = varargin{i}; 
     case {'topol','som_topol','sTopol'}, 
                        i=i+1; sTopol = varargin{i}; munits = prod(sTopol.msize); 
     case 'neigh',      i=i+1; neigh = varargin{i};
     case 'tracking',   i=i+1; tracking = varargin{i};
     case 'algorithm',  i=i+1; algorithm = varargin{i}; 
	 case 'distfun',	i=i+1; distfun = varargin{i};
     case 'init',       i=i+1; initalg = varargin{i};
     case 'training',   i=i+1; training = varargin{i}; 
	 case 'advanced',	i=i+1; advanced = varargin{i};
      % unambiguous values
     case {'hexa','rect'}, sTopol.lattice = varargin{i};
     case {'sheet','cyl','toroid'}, sTopol.shape = varargin{i}; 
     case {'gaussian','cutgauss','ep','bubble'}, neigh = varargin{i};
     case {'seq','batch','sompak'}, algorithm = varargin{i};
	 case {'sqEuclidean','correlation','ucorrelation'}, distfun = varargin{i};
     case {'small','normal','big'}, mapsize = varargin{i}; 
     case {'randinit','lininit'}, initalg = varargin{i};
     case {'short','default','long'}, training = varargin{i}; 
     otherwise argok=0; 
    end
  elseif isstruct(varargin{i}) && isfield(varargin{i},'type'), 
    switch varargin{i}(1).type, 
     case 'som_topol', sTopol = varargin{i};
     otherwise argok=0; 
    end
  else
    argok = 0; 
  end
  if ~argok, 
    disp(['(som_make) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the map struct

%% map size
%tole dodal Nejc Ilc
if isscalar(sTopol.msize)
    isRatio= (sTopol.msize < 0);
else
    isRatio=0;
end

if isempty(sTopol.msize) || ~prod(sTopol.msize) || isRatio, %zadnji clen dodal Nejc Ilc 
  if tracking>0, fprintf(1,'Determining map size...\n'); end
  
  %dodal Nejc Ilc
  %if sTopol.msize < 0, munits=0; end
  
  if ~munits,     
    sTemp = som_topol_struct('dlen',dlen);
    munits = prod(sTemp.msize);
    switch mapsize,
     case 'small', munits = max(9,ceil(munits/4));
     case 'big',   munits = munits*4;
     otherwise % nil
    end
  end
  %Spremenil Nejc Ilc - dodal za PCA
  sTemp = som_topol_struct('data',D,'munits',munits,'lattice',sTopol.lattice,'pca_val',pca_val);
  sTopol.msize = sTemp.msize;
  if tracking>0, 
    fprintf(1,' map size [%d, %d]\n',sTopol.msize(1), sTopol.msize(2));   
  end
end

% map struct - random positions of codebook vectors
sMap = som_map_struct(dim,sTopol,neigh,'mask',mask,'name',name, ...
                      'comp_names', comp_names, 'comp_norm', comp_norm); 

				  
				  
%% make sTrain struct from advanced settings, if exists
% sTrain struct - to modify advanced settings if neccessary
% If none of them provided, use defaults.
% Added by Nejc
sTrain = som_train_struct(sMap,'dlen',dlen,'algorithm',algorithm,'distfun',distfun);

if ~isempty(advanced)
	i=1; 
	while i<=length(advanced) 
		if ischar(advanced{i}) 
			switch advanced{i}
				case 'algorithm', i=i+1; algorithm = advanced{i}; sTrain=som_set(sTrain,'algorithm',algorithm);
				case 'initalg', i=i+1; initalg = advanced{i};
				case 'trainlen', i=i+1; training = advanced{i};
				case 'radius_ini', i=i+1; sTrain=som_set(sTrain,'radius_ini',advanced{i});
				case 'radius_fin', i=i+1; sTrain=som_set(sTrain,'radius_fin',advanced{i});
				case 'alpha_ini', i=i+1; sTrain=som_set(sTrain,'alpha_ini',advanced{i});
				case 'alpha_type', i=i+1; alpha_type=advanced{i}; sTrain=som_set(sTrain,'alpha_type',advanced{i});
				case 'sample_order', i=i+1; sample_order_type = advanced{i}; % random | ordered
				case 'neigh', i=i+1; neigh=advanced{i}; sTrain=som_set(sTrain,'neigh',neigh); sMap=som_set(sMap,'neigh',neigh);
				case 'dist_mode', i=i+1; dist_mode=advanced{i}; som_mode=strtok(dist_mode,'_');
				otherwise
					i=i+1;
			end
		end
		i=i+1;
	end
	if ~any(strcmp(som_mode,{'euc','dotprod','default'}))
		warning(['CreateSOM: Wrong dist_mode: ', dist_mode,', changing to default.']);
		som_mode = 'default';
	end
end

% added by Nejc Ilc
if any(strcmp(distfun,{'correlation','ucorrelation'}))
	if ~exist('som_batchtrain_corr','file')
		error('Function som_batchtrain_corr() is missing - cannot compute SOM for distance ''correlation''');
	end
	
	if strcmp(algorithm,'sompak')
		warning(['Algorithm ''',algorithm,''' not supported with distance ''',distfun,'''. Changing to batch mode!']);
		algorithm = 'batch';
	end
end
				  
				  
% function
if strcmp(algorithm,'sompak'), 
  algorithm = 'seq';
  func = 'sompak';
else
  func = algorithm;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization

if tracking>0, fprintf(1,'Initialization...\n'); end

switch initalg, 
	case 'randinit', sMap = som_randinit(D, sMap);
	case 'lininit', sMap = som_lininit(D, sMap, 'pca_vec',pca_vec,'pca_val',pca_val);
	%case 'lininit', sMap = som_nejcinit(D, sMap, 'pca_vec',pca_vec,'pca_val',pca_val);
end
sMap.trainhist(1) = som_set(sMap.trainhist(1),'data_name',data_name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% training

if tracking>0, fprintf(1,'Training using %s algorithm...\n',algorithm); end

% rough train
if tracking>0, fprintf(1,'Rough training phase...\n'); end


sTrain = som_train_struct(sTrain,'phase','rough');
sTrain = som_set(sTrain,'data_name',data_name,'distmode',dist_mode);


if isnumeric(training)
	sTrain.trainlen = training(1); 
else
  switch training, 
   case 'short', sTrain.trainlen = max(1,sTrain.trainlen/4);
   case 'long',  sTrain.trainlen = sTrain.trainlen*4;
  end
end
switch func,
 case 'seq',    
	switch distfun
		case {'sqEuclidean','default'}
			sMap = som_seqtrain(sMap,D,sTrain,'tracking',tracking,'mask',mask, 'sample_order', sample_order_type);
		case 'correlation'
			switch som_mode
				case {'euc','default'}
					sMap = som_seqtrain_corr(sMap,D,sTrain,'tracking',tracking,'mask',mask, 'sample_order', sample_order_type);
				case 'dotprod'
					sMap = som_seqtrain_dotprodc(sMap,D,sTrain,'tracking',tracking,'mask',mask, 'sample_order', sample_order_type);
			end
		case 'ucorrelation'
			switch som_mode
				case {'euc','default'}
					sMap = som_seqtrain_ucorr(sMap,D,sTrain,'tracking',tracking,'mask',mask, 'sample_order', sample_order_type);
				case 'dotprod'
					sMap = som_seqtrain_dotprod(sMap,D,sTrain,'tracking',tracking,'mask',mask, 'sample_order', sample_order_type);
			end
		otherwise
			error(['Wrong distance: ',distfun]);
	end
 case 'sompak', sMap = som_sompaktrain(sMap,D,sTrain,'tracking',tracking,'mask',mask);
 case 'batch',
	switch distfun
		case {'sqEuclidean','default'}
			sMap = som_batchtrain(sMap,D,sTrain,'tracking',tracking,'mask',mask);
		case 'correlation'
			switch som_mode
				case {'euc','default'}
					sMap = som_batchtrain_corr(sMap,D,sTrain,'tracking',tracking,'mask',mask); 
				case 'dotprod'
					sMap = som_batchtrain_dotprodc(sMap,D,sTrain,'tracking',tracking,'mask',mask);
			end
		case 'ucorrelation'
			switch som_mode
				case {'euc','default'}
					sMap = som_batchtrain_ucorr(sMap,D,sTrain,'tracking',tracking,'mask',mask);
				case 'dotprod'
					sMap = som_batchtrain_dotprod(sMap,D,sTrain,'tracking',tracking,'mask',mask);
			end
		otherwise
			error(['Wrong distance: ',distfun]);
	end
otherwise
	error(['Wrong algorithm: ',func]);
end


%---------------------------------------------------------------
% finetune
if tracking>0, fprintf(1,'Finetuning phase...\n'); end
sTrain_tmp = som_train_struct(sMap,sTrain,'dlen',dlen,'phase','finetune');
sTrain = som_set(sTrain_tmp, 'radius_ini',sTrain.radius_fin,'radius_fin',1);
%som_set(sTrain,'data_name',data_name,'algorithm',algorithm,'alpha_type',alpha_type,'distfun',distfun);

if isnumeric(training), sTrain.trainlen = training(2); 
else
  switch training, 
   case 'short', sTrain.trainlen = max(1,sTrain.trainlen/4);
   case 'long',  sTrain.trainlen = sTrain.trainlen*4;
  end
end

switch func,
 case 'seq',    
	switch distfun
		case {'sqEuclidean','default'}
			sMap = som_seqtrain(sMap,D,sTrain,'tracking',tracking,'mask',mask, 'sample_order', sample_order_type);
			
		case 'correlation'
			switch som_mode
				case {'euc','default'}
					sMap = som_seqtrain_corr(sMap,D,sTrain,'tracking',tracking,'mask',mask, 'sample_order', sample_order_type);
				case 'dotprod'
					sMap = som_seqtrain_dotprodc(sMap,D,sTrain,'tracking',tracking,'mask',mask, 'sample_order', sample_order_type);
			end
		case 'ucorrelation'
			switch som_mode
				case {'euc','default'}
					sMap = som_seqtrain_ucorr(sMap,D,sTrain,'tracking',tracking,'mask',mask, 'sample_order', sample_order_type);
				case 'dotprod'
					sMap = som_seqtrain_dotprod(sMap,D,sTrain,'tracking',tracking,'mask',mask, 'sample_order', sample_order_type);
			end
		otherwise
			error(['Wrong distance: ',distfun]);
	end
 case 'sompak', sMap = som_sompaktrain(sMap,D,sTrain,'tracking',tracking,'mask',mask);
 case 'batch',
	switch distfun
		case {'sqEuclidean','default'}
			sMap = som_batchtrain(sMap,D,sTrain,'tracking',tracking,'mask',mask);
			
		case 'correlation'
			switch som_mode
				case {'euc','default'}
					sMap = som_batchtrain_corr(sMap,D,sTrain,'tracking',tracking,'mask',mask); 
				case 'dotprod'
					sMap = som_batchtrain_dotprodc(sMap,D,sTrain,'tracking',tracking,'mask',mask);
			end
		case 'ucorrelation'
			switch som_mode
				case {'euc','default'}
					sMap = som_batchtrain_ucorr(sMap,D,sTrain,'tracking',tracking,'mask',mask);
				case 'dotprod'
					sMap = som_batchtrain_dotprod(sMap,D,sTrain,'tracking',tracking,'mask',mask);
			end
		otherwise
			error(['Wrong distance: ',distfun]);
	end
otherwise
	error(['Wrong algorithm: ',func]);
end
%---------------------------------------------------------------


% quality
if tracking>0, 
  [mqe,tge] = som_quality(sMap,D);
  fprintf(1,'Final quantization error: %5.3f\n',mqe)
  fprintf(1,'Final topographic error:  %5.3f\n',tge)
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

