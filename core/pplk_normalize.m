function sD = pplk_normalize(sD,method,comps)
%based on: SOM_NORMALIZE (Re)normalize data or add new normalizations.
%
% sS = som_normalize(sS,[method],[comps])
%
%   sS = som_normalize(sD)
%   sS = som_normalize(sS,sNorm)
%    D = som_normalize(D,'var')
%   sS = som_normalize(sS,'histC',[1:3 10])
%
%  Input and output arguments ([]'s are optional):
%   sS                The data to which the normalization is applied.
%                     The modified and updated data is returned.
%            (struct) data or map struct
%            (matrix) data matrix (a matrix is also returned)
%   [method]          The normalization method(s) to add/use. If missing,
%                     or an empty variable ('') is given, the
%                     normalizations in sS are used.
%            (string) identifier for a normalization method to be added:
%                     'var', 
%                     'range', 
%                     'log', 
%                     'logistic', 
%                     'histD', 
%                     'histC',
%                     'norm',
%                     'zscore',
%                     'propor'.
%
%            (struct) Normalization struct, or an array of such.
%                     Alternatively, a map/data struct can be given
%                     in which case its '.comp_norm' field is used
%                     (see below).
%            (cell array) Of normalization structs. Typically, the
%                     '.comp_norm' field of a map/data struct. The
%                     length of the array must be equal to data dimension.
%            (cellstr array) norm and denorm operations in a cellstr array
%                     which are evaluated with EVAL command with variable
%                     name 'x' reserved for the variable.
%   [comps]  (vector) the components to which the normalization is
%                     applied, default is [1:dim] ie. all components
%
% For more help, try 'type som_normalize' or check out online documentation.
% See also SOM_DENORMALIZE, SOM_NORM_VARIABLE, SOM_INFO.

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% som_normalize
%
% PURPOSE
%
% Add/apply/redo normalization on data structs/sets.
%
% SYNTAX
%
%  sS = som_normalize(sS)
%  sS = som_normalize(sS,method)
%   D = som_normalize(D,sNorm)
%  sS = som_normalize(sS,csNorm)
%  sS = som_normalize(...,comps)
%
% DESCRIPTION
%
% This function is used to (initialize and) add, redo and apply
% normalizations on data/map structs/sets. If a data/map struct is given,
% the specified normalizations are added to the '.comp_norm' field of the
% struct after ensuring that all normalizations specified therein have
% status 'done'. SOM_NORMALIZE actually uses function SOM_NORM_VARIABLE
% to handle the normalization operations, and only handles the data
% struct/set specific stuff itself.
%
% The different normalization methods are listed below. For more
% detailed descriptions, see SOM_NORM_VARIABLE.
%
%   method     description
%   'none'     skip noramlization, data unchanged
%   'var'      Variance is normalized to one (linear operation).
%   'range'    Values are normalized between [0,1] (linear operation).
%   'log'      Natural logarithm is applied to the values:
%                xnew = log(x-m+1)
%              where m = min(x).
%   'logistic' Logistic or softmax trasformation which scales all
%              possible values between [0,1].
%   'histD'    Histogram equalization, values scaled between [0,1].
%   'histC'    Approximate histogram equalization with partially
%              linear operations. Values scaled between [0,1].
%   'eval'     freeform operations
%   'zscore'   same as 'var'! substract data vector the data mean and divide  
%              by stdev of data (alse called zero-mean-unit-variance)
%   'norm'     each data vector becomes of unit length (||data_i|| = 1)
%	'propor'   proportional scaling to [0,1].
%			   Same as 'range' except that it preserves ratios between
%			   variables.
%	'maxdist'  normalize on max pair-wise distance.
%
% To enable undoing and applying the exactly same normalization to
% other data sets, normalization information is saved into a
% normalization struct, which has the fields:
%
%   .type   ; struct type, ='som_norm'
%   .method ; normalization method, a string
%   .params ; normalization parameters
%   .status ; string: 'uninit', 'undone' or 'done'
%
% Normalizations are always one-variable operations. In the data and map
% structs the normalization information for each component is saved in the
% '.comp_norm' field, which is a cell array of length dim. Each cell
% contains normalizations for one vector component in a struct array of
% normalization structs. Each component may have different amounts of
% different kinds of normalizations. Typically, all normalizations are
% either 'undone' or 'done', but in special situations this may not be the
% case. The easiest way to check out the status of the normalizations is to
% use function SOM_INFO, e.g. som_info(sS,3)
%
% REQUIRED INPUT ARGUMENTS
%
%   sS                The data to which the normalization is applied.
%            (struct) Data or map struct. Before adding any new
%                     normalizations, it is ensured that the
%                     normalizations for the specified components in the
%                     '.comp_norm' field have status 'done'.
%            (matrix) data matrix
%
% OPTIONAL INPUT ARGUMENTS
%
%   method            The normalization(s) to add/use. If missing,
%                     or an empty variable ('' or []) is given, the
%                     normalizations in the data struct are used.
%            (string) Identifier for a normalization method to be added:
%                     'var', 'range', 'log', 'logistic', 'histD' or 'histC'. The
%                     same method is applied to all specified components
%                     (given in comps). The normalizations are first
%                     initialized (for each component separately, of
%                     course) and then applied.
%            (struct) Normalization struct, or an array of structs, which
%                     is applied to all specified components. If the
%                     '.status' field of the struct(s) is 'uninit',
%                     the normalization(s) is initialized first.
%                     Alternatively, the struct may be map or data struct
%                     in which case its '.comp_norm' field is used
%                     (see the cell array option below).
%            (cell array) In practice, the '.comp_norm' field of
%                     a data/map struct. The length of the array
%                     must be equal to the dimension of the given
%                     data set (sS). Each cell contains the
%                     normalization(s) for one component. Only the
%                     normalizations listed in comps argument are
%                     applied though.
%            (cellstr array) norm and denorm operations in a cellstr array
%                     which are evaluated with EVAL command with variable
%                     name 'x' reserved for the variable.
%
%   comps    (vector) The components to which the normalization(s) is
%                     applied. Default is to apply to all components.
%
% OUTPUT ARGUMENTS
%
%   sS                Modified and/or updated data.
%            (struct) If a struct was given as input argument, the
%                     same struct is returned with normalized data and
%                     updated '.comp_norm' fields.
%            (matrix) If a matrix was given as input argument, the
%                     normalized data matrix is returned.
%
% EXAMPLES
%
%  To add (initialize and apply) a normalization to a data struct:
%
%    sS = som_normalize(sS,'var');
%
%  This uses 'var'-method to all components. To add a method only to
%  a few selected components, use the comps argument:
%
%    sS = som_normalize(sS,'log',[1 3:5]);
%
%  To ensure that all normalization operations have indeed been done:
%
%    sS = som_normalize(sS);
%
%  The same for only a few components:
%
%    sS = som_normalize(sS,'',[1 3:5]);
%
%  To apply the normalizations of a data struct sS to a new data set D:
%
%    D = som_normalize(D,sS);
%  or
%    D = som_normalize(D,sS.comp_norm);
%
%  To normalize a data set:
%
%    D = som_normalize(D,'histD');
%
%  Note that in this case the normalization information is lost.
%
%  To check out the status of normalization in a struct use SOM_INFO:
%
%    som_info(sS,3)
%
%
% SEE ALSO
%
%  som_denormalize    Undo normalizations of a data struct/set.
%  som_norm_variable  Normalization operations for a set of scalar values.
%  som_info           User-friendly information of SOM Toolbox structs.

% Copyright (c) 1998-2000 by the SOM toolbox programming team.
% http://www.cis.hut.fi/projects/somtoolbox/

% Version 2.0beta juuso 151199 150500
% Modified by Nejc Ilc - added 'norm', 'zscore', 'propor', 'maxdist'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

callDir=chdir(pplk_homeDir());
SOMTOOLBOX_PATH = ['..',filesep,'libs',filesep,'somtoolbox_pplk'];
addpath(SOMTOOLBOX_PATH);

%% check arguments

narginchk(1, 3);  % check no. of input arguments is correct

% sD
struct_mode = isstruct(sD);
if struct_mode,
    switch sD.type
        case 'som_map', D = sD.codebook;
        case 'som_data', D = sD.data;
        otherwise, error('Illegal struct.')
    end
else
    D = sD;
end
[dlen dim] = size(D);

% comps
if nargin<3 || (ischar(comps) && strcmp(comps,'all')),
    comps = [1:dim];
end
if isempty(comps), chdir(callDir);return; end
if size(comps,1)>1, comps = comps'; end  % make it a row vector

% method
csNorm = cell(dim,1);
if nargin<2 || isempty(method)
    if ~struct_mode,
        warning('No normalization method given. Data left unchanged.');
        chdir(callDir);
        return;
    end
    method = '';
else
    % check out the given method
    % (and if necessary, copy it for each specified component)
    if strcmpi(method,'none')
        chdir(callDir);
        return;
    end
    
    if ischar(method),
        switch method,
            case {'var','range','log','histD','histC','logistic','norm','zscore','propor','maxdist'},
                sN = som_set('som_norm','method',method);
            otherwise,
                error(['Unrecognized method: ' method]);
        end
        for i=comps, csNorm{i} = sN; end
    elseif isstruct(method),
        switch method(1).type,
            case {'som_map','som_data'}, csNorm = method(1).comp_norm;
            case {'som_norm'}, for i=comps, csNorm{i} = method; end
            otherwise,
                error('Invalid struct given as normalization method.')
        end
    elseif iscellstr(method),
        [dummy,sN] = som_norm_variable(1,method,'init');
        for i=comps, csNorm{i} = sN; end
    elseif iscell(method),
        csNorm = method;
    else
        error('Illegal method argument.')
    end
    % check the size of csNorm is the same as data dimension
    if length(csNorm) ~= dim,
        error('Given number of normalizations does not match data dimension.')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialize

% make sure all the current normalizations for current
% components have been done
if struct_mode,
    alldone = 1;
    for i = comps,
        for j=1:length(sD.comp_norm{i}),
            sN = sD.comp_norm{i}(j);
            if ~strcmp(sN.status,'done'),
                alldone = 0;
                [x,sN] = som_norm_variable(D(:,i), sN, 'do');
                D(:,i) = x;
                sD.comp_norm{i}(j) = sN;
            end
        end
    end
    if isempty(method),
        if alldone,
            warning('No ''undone'' normalizations found. Data left unchanged.');
        else
            fprintf(1,'Normalizations have been redone.\n');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% action

if strcmpi(method,'maxdist')
    % we take all the variables into account, so it is easier to consider
    % it seperately.
    p=norm_maxdist_init(D);
    D=norm_maxdist_do(D,p);

elseif strcmpi(method,'propor')
    % we take all the variables into account, so it is easier to consider
    % it seperately.
    p=norm_propor_init(D);
    D=norm_propor_do(D,p);
    
elseif strcmpi(method,'norm')
    p=norm_norm_init(D);
    D=norm_norm_do(D,p);
    
elseif strcmpi(method,'zscore')
    p=norm_zscore_init(D);
    D=norm_zscore_do(D,p);
    
else
    
    % add the new normalizations to the old ones
    for i = comps,
        if ~isempty(csNorm{i}),
            [x,sN] = som_norm_variable(D(:,i), csNorm{i}, 'do');
            D(:,i) = x;
            if struct_mode,
                if isempty(sD.comp_norm{i}), sD.comp_norm{i} = sN;
                else sD.comp_norm{i} = [sD.comp_norm{i}, sN]; end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% output

if struct_mode,
    switch sD.type
        case 'som_map', sD.codebook = D;
        case 'som_data', sD.data = D;
        otherwise, error('Illegal struct.')
    end
else
    sD = D;
end

%remove SOM Toolbox from MATLAB path
rmpath(SOMTOOLBOX_PATH);

chdir(callDir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,sNorm] = som_norm_variable(x, method, operation)

%SOM_NORM_VARIABLE Normalize or denormalize a scalar variable.
%
% [x,sNorm] = som_norm_variable(x, method, operation)
%
%   xnew = som_norm_variable(x,'var','do');
%   [dummy,sN] = som_norm_variable(x,'log','init');
%   [xnew,sN]  = som_norm_variable(x,sN,'do');
%   xorig      = som_norm_variable(xnew,sN,'undo');
%
%  Input and output arguments:
%   x         (vector) a set of values of a scalar variable for
%                      which the (de)normalization is performed.
%                      The processed values are returned.
%   method    (string) identifier for a normalization method: 'var',
%                      'range', 'log', 'logistic', 'histD', or 'histC'.
%                      A normalization struct with default values is created.
%             (struct) normalization struct, or an array of such
%             (cellstr) first string gives normalization operation, and the
%                      second gives denormalization operation, with x
%                      representing the variable, for example:
%                      {'x+2','x-2}, or {'exp(-x)','-log(x)'} or {'round(x)'}.
%                      Note that in the last case, no denorm operation is
%                      defined.
%   operation (string) the operation to be performed: 'init', 'do' or 'undo'
%
%   sNorm     (struct) updated normalization struct/struct array
%
% For more help, try 'type som_norm_variable' or check out online documentation.
% See also SOM_NORMALIZE, SOM_DENORMALIZE.
%
%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% som_norm_variable
%
% PURPOSE
%
% Initialize, apply and undo normalizations on a given vector of
% scalar values.
%
% SYNTAX
%
%  xnew = som_norm_variable(x,method,operation)
%  xnew = som_norm_variable(x,sNorm,operation)
%  [xnew,sNorm] = som_norm_variable(...)
%
% DESCRIPTION
%
% This function is used to initialize, apply and undo normalizations
% on scalar variables. It is the low-level function that upper-level
% functions SOM_NORMALIZE and SOM_DENORMALIZE utilize to actually (un)do
% the normalizations.
%
% Normalizations are typically performed to control the variance of
% vector components. If some vector components have variance which is
% significantly higher than the variance of other components, those
% components will dominate the map organization. Normalization of
% the variance of vector components (method 'var') is used to prevent
% that. In addition to variance normalization, other methods have
% been implemented as well (see list below).
%
% Usually normalizations convert the variable values so that they no
% longer make any sense: the values are still ordered, but their range
% may have changed so radically that interpreting the numbers in the
% original context is very hard. For this reason all implemented methods
% are (more or less) revertible. The normalizations are monotonic
% and information is saved so that they can be undone. Also, the saved
% information makes it possible to apply the EXACTLY SAME normalization
% to another set of values. The normalization information is determined
% with 'init' operation, while 'do' and 'undo' operations are used to
% apply or revert the normalization.
%
% The normalization information is saved in a normalization struct,
% which is returned as the second argument of this function. Note that
% normalization operations may be stacked. In this case, normalization
% structs are positioned in a struct array. When applied, the array is
% gone through from start to end, and when undone, in reverse order.
%
%    method  description
%
%    'var'   Variance normalization. A linear transformation which
%            scales the values such that their variance=1. This is
%            convenient way to use Mahalanobis distance measure without
%            actually changing the distance calculation procedure.
%
%    'range' Normalization of range of values. A linear transformation
%            which scales the values between [0,1].
%
%    'log'   Logarithmic normalization. In many cases the values of
%            a vector component are exponentially distributed. This
%            normalization is a good way to get more resolution to
%            (the low end of) that vector component. What this
%            actually does is a non-linear transformation:
%               x_new = log(x_old - m + 1)
%            where m=min(x_old) and log is the natural logarithm.
%            Applying the transformation to a value which is lower
%            than m-1 will give problems, as the result is then complex.
%            If the minimum for values is known a priori,
%            it might be a good idea to initialize the normalization with
%              [dummy,sN] = som_norm_variable(minimum,'log','init');
%            and normalize only after this:
%              x_new = som_norm_variable(x,sN,'do');
%
%    'logistic' or softmax normalization. This normalization ensures
%            that all values in the future, too, are within the range
%            [0,1]. The transformation is more-or-less linear in the
%            middle range (around mean value), and has a smooth
%            nonlinearity at both ends which ensures that all values
%            are within the range. The data is first scaled as in
%            variance normalization:
%               x_scaled = (x_old - mean(x_old))/std(x_old)
%            and then transformed with the logistic function
%               x_new = 1/(1+exp(-x_scaled))
%
%    'histD' Discrete histogram equalization. Non-linear. Orders the
%            values and replaces each value by its ordinal number.
%            Finally, scales the values such that they are between [0,1].
%            Useful for both discrete and continuous variables, but as
%            the saved normalization information consists of all
%            unique values of the initialization data set, it may use
%            considerable amounts of memory. If the variable can get
%            more than a few values (say, 20), it might be better to
%            use 'histC' method below. Another important note is that
%            this method is not exactly revertible if it is applied
%            to values which are not part of the original value set.
%
%    'histC' Continuous histogram equalization. Actually, a partially
%            linear transformation which tries to do something like
%            histogram equalization. The value range is divided to
%            a number of bins such that the number of values in each
%            bin is (almost) the same. The values are transformed
%            linearly in each bin. For example, values in bin number 3
%            are scaled between [3,4[. Finally, all values are scaled
%            between [0,1]. The number of bins is the square root
%            of the number of unique values in the initialization set,
%            rounded up. The resulting histogram equalization is not
%            as good as the one that 'histD' makes, but the benefit
%            is that it is exactly revertible - even outside the
%            original value range (although the results may be funny).
%
%    'eval'  With this method, freeform normalization operations can be
%            specified. The parameter field contains strings to be
%            evaluated with 'eval' function, with variable name 'x'
%            representing the variable itself. The first string is
%            the normalization operation, and the second is a
%            denormalization operation. If the denormalization operation
%            is empty, it is ignored.
%
% INPUT ARGUMENTS
%
%   x          (vector) The scalar values to which the normalization
%                       operation is applied.
%
%   method              The normalization specification.
%              (string) Identifier for a normalization method: 'var',
%                       'range', 'log', 'logistic', 'histD' or 'histC'.
%                       Corresponding default normalization struct is created.
%              (struct) normalization struct
%              (struct array) of normalization structs, applied to
%                       x one after the other
%              (cellstr) of length
%              (cellstr array) first string gives normalization operation, and
%                       the second gives denormalization operation, with x
%                       representing the variable, for example:
%                       {'x+2','x-2}, or {'exp(-x)','-log(x)'} or {'round(x)'}.
%                       Note that in the last case, no denorm operation is
%                       defined.
%
%               note: if the method is given as struct(s), it is
%                     applied (done or undone, as specified by operation)
%                     regardless of what the value of '.status' field
%                     is in the struct(s). Only if the status is
%                     'uninit', the undoing operation is halted.
%                     Anyhow, the '.status' fields in the returned
%                     normalization struct(s) is set to approriate value.
%
%   operation  (string) The operation to perform: 'init' to initialize
%                       the normalization struct, 'do' to perform the
%                       normalization, 'undo' to undo the normalization,
%                       if possible. If operation 'do' is given, but the
%                       normalization struct has not yet been initialized,
%                       it is initialized using the given data (x).
%
% OUTPUT ARGUMENTS
%
%   x        (vector) Appropriately processed values.
%
%   sNorm    (struct) Updated normalization struct/struct array. If any,
%                     the '.status' and '.params' fields are updated.
%
% EXAMPLES
%
%  To initialize and apply a normalization on a set of scalar values:
%
%    [x_new,sN] = som_norm_variable(x_old,'var','do');
%
%  To just initialize, use:
%
%    [dummy,sN] = som_norm_variable(x_old,'var','init');
%
%  To undo the normalization(s):
%
%    x_orig = som_norm_variable(x_new,sN,'undo');
%
%  Typically, normalizations of data structs/sets are handled using
%  functions SOM_NORMALIZE and SOM_DENORMALIZE. However, when only the
%  values of a single variable are of interest, SOM_NORM_VARIABLE may
%  be useful. For example, assume one wants to apply the normalization
%  done on a component (i) of a data struct (sD) to a new set of values
%  (x) of that component. With SOM_NORM_VARIABLE this can be done with:
%
%    x_new = som_norm_variable(x,sD.comp_norm{i},'do');
%
%  Now, as the normalizations in sD.comp_norm{i} have already been
%  initialized with the original data set (presumably sD.data),
%  the EXACTLY SAME normalization(s) can be applied to the new values.
%  The same thing can be done with SOM_NORMALIZE function, too:
%
%    x_new = som_normalize(x,sD.comp_norm{i});
%
%  Or, if the new data set were in variable D - a matrix of same
%  dimension as the original data set:
%
%    D_new = som_normalize(D,sD,i);
%
% SEE ALSO
%
%  som_normalize    Add/apply/redo normalizations for a data struct/set.
%  som_denormalize  Undo normalizations of a data struct/set.

% Copyright (c) 1998-2000 by the SOM toolbox programming team.
% http://www.cis.hut.fi/projects/somtoolbox/

% Version 2.0beta juuso 151199 170400 150500

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check arguments

narginchk(3, 3);  % check no. of input arguments is correct

% method
sNorm = [];
if ischar(method)
    if any(strcmp(method,{'var','range','log','logistic','histD','histC'})),
        sNorm = som_set('som_norm','method',method);
    else
        method = cellstr(method);
    end
end
if iscell(method),
    if length(method)==1 & isstruct(method{1}), sNorm = method{1};
    else
        if length(method)==1 | isempty(method{2}), method{2} = 'x'; end
        sNorm = som_set('som_norm','method','eval','params',method);
    end
else
    sNorm = method;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% action

order = [1:length(sNorm)];
if length(order)>1 & strcmp(operation,'undo'), order = order(end:-1:1); end

for i=order,
    
    % initialize
    if strcmp(operation,'init') | ...
            (strcmp(operation,'do') & strcmp(sNorm(i).status,'uninit')),
        
        % case method = 'hist'
        if strcmp(sNorm(i).method,'hist'),
            inds = find(~isnan(x) & ~isinf(x));
            if length(unique(x(inds)))>20, sNorm(i).method = 'histC';
            else sNorm{i}.method = 'histD'; end
        end
        
        switch(sNorm(i).method),
            case 'var',   params = norm_variance_init(x);
            case 'range', params = norm_scale01_init(x);
            case 'log',   params = norm_log_init(x);
            case 'logistic', params = norm_logistic_init(x);
            case 'histD', params = norm_histeqD_init(x);
            case 'histC', params = norm_histeqC_init(x);
            case 'eval',  params = sNorm(i).params;
            case 'norm', params = norm_norm_init(x);
            case 'zscore', params = norm_zscore_init(x);
            case 'propor', params = norm_propor_init(x);
            otherwise,
                error(['Unrecognized method: ' sNorm(i).method]);
        end
        sNorm(i).params = params;
        sNorm(i).status = 'undone';
    end
    
    % do / undo
    if strcmp(operation,'do'),
        switch(sNorm(i).method),
            case 'var',   x = norm_scale_do(x,sNorm(i).params);
            case 'range', x = norm_scale_do(x,sNorm(i).params);
            case 'log',   x = norm_log_do(x,sNorm(i).params);
            case 'logistic', x = norm_logistic_do(x,sNorm(i).params);
            case 'histD', x = norm_histeqD_do(x,sNorm(i).params);
            case 'histC', x = norm_histeqC_do(x,sNorm(i).params);
            case 'eval',  x = norm_eval_do(x,sNorm(i).params);
            case 'norm', params = norm_norm_do(x,sNorm(i).params);
            case 'zscore', params = norm_zscore_do(x,sNorm(i).params);
            case 'propor', params = norm_propor_do(x,sNorm(i).params);
            otherwise,
                error(['Unrecognized method: ' sNorm(i).method]);
        end
        sNorm(i).status = 'done';
        
    elseif strcmp(operation,'undo'),
        
        if strcmp(sNorm(i).status,'uninit'),
            warning('Could not undo: uninitialized normalization struct.')
            break;
        end
        switch(sNorm(i).method),
            case 'var',   x = norm_scale_undo(x,sNorm(i).params);
            case 'range', x = norm_scale_undo(x,sNorm(i).params);
            case 'log',   x = norm_log_undo(x,sNorm(i).params);
            case 'logistic', x = norm_logistic_undo(x,sNorm(i).params);
            case 'histD', x = norm_histeqD_undo(x,sNorm(i).params);
            case 'histC', x = norm_histeqC_undo(x,sNorm(i).params);
            case 'eval',  x = norm_eval_undo(x,sNorm(i).params);
            case 'norm', params = norm_norm_undo(x,sNorm(i).params);
            case 'zscore', params = norm_zscore_undo(x,sNorm(i).params);
            case 'propor', params = norm_propor_undo(x,sNorm(i).params);
            otherwise,
                error(['Unrecognized method: ' sNorm(i).method]);
        end
        sNorm(i).status = 'undone';
        
    elseif ~strcmp(operation,'init'),
        
        error(['Unrecognized operation: ' operation])
        
    end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunctions

% linear scaling

function p = norm_variance_init(x)
inds = find(~isnan(x) & isfinite(x));
p = [mean(x(inds)), std(x(inds))];
if p(2) == 0, p(2) = 1; end
%end of norm_variance_init

function p = norm_scale01_init(x)
inds = find(~isnan(x) & isfinite(x));
mi = min(x(inds));
ma = max(x(inds));
if mi == ma, p = [mi, 1]; else p = [mi, ma-mi]; end
%end of norm_scale01_init

function x = norm_scale_do(x,p)
x = (x - p(1)) / p(2);
% end of norm_scale_do

function x = norm_scale_undo(x,p)
x = x * p(2) + p(1);
% end of norm_scale_undo

% logarithm

function p = norm_log_init(x)
inds = find(~isnan(x) & isfinite(x));
p = min(x(inds));
% end of norm_log_init

function x = norm_log_do(x,p)
x = log(x - p +1);
% if any(~isreal(x)), ok = 0; end
% end of norm_log_do

function x = norm_log_undo(x,p)
x = exp(x) -1 + p;
% end of norm_log_undo

% logistic

function p = norm_logistic_init(x)
inds = find(~isnan(x) & isfinite(x));
p = [mean(x(inds)), std(x(inds))];
if p(2)==0, p(2) = 1; end
% end of norm_logistic_init

function x = norm_logistic_do(x,p)
x = (x-p(1))/p(2);
x = 1./(1+exp(-x));
% end of norm_logistic_do

function x = norm_logistic_undo(x,p)
x = log(x./(1-x));
x = x*p(2)+p(1);
% end of norm_logistic_undo

% histogram equalization for discrete values

function p = norm_histeqD_init(x)
inds = find(~isnan(x) & ~isinf(x));
p = unique(x(inds));
% end of norm_histeqD_init

function x = norm_histeqD_do(x,p)
bins = length(p);
inds = find(~isnan(x) & ~isinf(x))';
for i = inds,
    [dummy ind] = min(abs(x(i) - p));
    % data item closer to the left-hand bin wall is indexed after RH wall
    if x(i) > p(ind) & ind < bins,
        x(i) = ind + 1;
    else
        x(i) = ind;
    end
end
x = (x-1)/(bins-1); % normalization between [0,1]
% end of norm_histeqD_do

function x = norm_histeqD_undo(x,p)
bins = length(p);
x = round(x*(bins-1)+1);
inds = find(~isnan(x) & ~isinf(x));
x(inds) = p(x(inds));
% end of norm_histeqD_undo

% histogram equalization with partially linear functions

function p = norm_histeqC_init(x)
% investigate x
inds = find(~isnan(x) & ~isinf(x));
samples = length(inds);
xs = unique(x(inds));
mi = xs(1);
ma = xs(end);
% decide number of limits
lims = ceil(sqrt(length(xs))); % 2->2,100->10,1000->32,10000->100
% decide limits
if lims==1,
    p = [mi, mi+1];
    lims = 2;
elseif lims==2,
    p = [mi, ma];
else
    p = zeros(lims,1);
    p(1) = mi;
    p(end) = ma;
    binsize = zeros(lims-1,1); b = 1; avebinsize = samples/(lims-1);
    for i=1:(length(xs)-1),
        binsize(b) = binsize(b) + sum(x==xs(i));
        if binsize(b) >= avebinsize,
            b = b + 1;
            p(b) = (xs(i)+xs(i+1))/2;
        end
        if b==(lims-1),
            binsize(b) = samples-sum(binsize); break;
        else
            avebinsize = (samples-sum(binsize))/(lims-1-b);
        end
    end
end
% end of norm_histeqC_init

function x = norm_histeqC_do(x,p)
xnew = x;
lims = length(p);
% handle values below minimum
r = p(2)-p(1);
inds = find(x<=p(1) & isfinite(x));
if any(inds), xnew(inds) = 0-(p(1)-x(inds))/r; end
% handle values above maximum
r = p(end)-p(end-1);
inds = find(x>p(end) & isfinite(x));
if any(inds), xnew(inds) = lims-1+(x(inds)-p(end))/r; end
% handle all other values
for i=1:(lims-1),
    r0 = p(i); r1 = p(i+1); r = r1-r0;
    inds = find(x>r0 & x<=r1);
    if any(inds), xnew(inds) = i-1+(x(inds)-r0)/r; end
end
% scale so that minimum and maximum correspond to 0 and 1
x = xnew/(lims-1);
% end of norm_histeqC_do

function x = norm_histeqC_undo(x,p)
xnew = x;
lims = length(p);
% scale so that 0 and 1 correspond to minimum and maximum
x = x*(lims-1);

% handle values below minimum
r = p(2)-p(1);
inds = find(x<=0 & isfinite(x));
if any(inds), xnew(inds) = x(inds)*r + p(1); end
% handle values above maximum
r = p(end)-p(end-1);
inds = find(x>lims-1 & isfinite(x));
if any(inds), xnew(inds) = (x(inds)-(lims-1))*r+p(end); end
% handle all other values
for i=1:(lims-1),
    r0 = p(i); r1 = p(i+1); r = r1-r0;
    inds = find(x>i-1 & x<=i);
    if any(inds), xnew(inds) = (x(inds)-(i-1))*r + r0; end
end
x = xnew;
% end of norm_histeqC_undo

% eval

function p = norm_eval_init(method)
p = method;
%end of norm_eval_init

function x = norm_eval_do(x,p)
x_tmp = eval(p{1});
if size(x_tmp,1)>=1 & size(x,1)>=1 & ...
        size(x_tmp,2)==1 & size(x,2)==1,
    x = x_tmp;
end
%end of norm_eval_do

function x = norm_eval_undo(x,p)
x_tmp = eval(p{2});
if size(x_tmp,1)>=1 & size(x,1)>=1 & ...
        size(x_tmp,2)==1 & size(x,2)==1,
    x = x_tmp;
end
%end of norm_eval_undo

% added by Nejc

% Normalization of data vectors to have unit length, i.e., ||data_i|| = 1
function p = norm_norm_init(D)
% compute norms of each data vector

p = sqrt(sum(D.^2,2));
%end of norm_norm_init

function D = norm_norm_do(D,p)
D = bsxfun(@rdivide,D, p);
% end of norm_norm_do

function D = norm_norm_undo(D,p)
D = bsxfun(@times,D, p);
% end of norm_norm_undo


% Z-score
function p = norm_zscore_init(D)
% compute mean and std of data
p = [mean(D,1); std(D,0,1)];
%end of norm_norm_init

function D = norm_zscore_do(D,p)
tmp = bsxfun(@minus,D,p(1,:));
D = bsxfun(@rdivide,tmp,p(2,:));
% end of norm_norm_do

function D = norm_zscore_undo(D,p)
D = bsxfun(@times,D,p(2,:));
D = bsxfun(@plus,D,p(1,:));
% end of norm_norm_undo


% Normalization on interval [0,1] perserving ratios between variables (unlike 'range').
% After normalization, variable with the highest dynamic range (max-min)
% will still have the highest max-min value
% (e.g., rectangular (2 x 1) becomes rectangular (1 x 0.5) and not a square (1 x 1))

function p = norm_propor_init(D)
%find the largest max-min and compute scaling coefficient
[n,d]=size(D);
inds=zeros(n,1);

ma=zeros(1,d);
mi=zeros(1,2);

for dim=1:d
    inds = find(~isnan(D(:,dim)) & isfinite(D(:,dim)));
    ma(dim)=max(D(inds,dim));
    mi(dim)=min(D(inds,dim));
end
% max difference between max and min of particular variable
r = max(ma-mi);

p = {mi,r};

%end of norm_propor_init

function D = norm_propor_do(D,p)
D = (D - repmat(p{1},size(D,1),1)) / p{2};
% end of norm_propor_do

function D = norm_propor_undo(D,p)
D = D * p(2) + repmat(p(1),size(D),1);
% end of norm_propor_undo


% Normalization on maximum pair-wise distance between points.
function p = norm_maxdist_init(D)
DIST = squareform(pdist(D,'euclidean'));
maxDi_val = max(DIST,[],1);
[maxD,maxD_ind] = max(maxDi_val);
d1 = maxD_ind;
p = [maxD,maxD_ind];


function D = norm_maxdist_do(D,p)
maxD = p(1);
d1 = p(2);
% move all data, so that d1 is in the origin (0,0)
dataNorm = bsxfun(@minus,D,D(d1,:));
% divide all data points by maxD
dataNorm = bsxfun(@rdivide, dataNorm, maxD);
% move back to original position of d1
D = bsxfun(@plus,dataNorm,D(d1,:));


function D = norm_maxdist_undo(D,p)
maxD = p(1);
d1 = p(2);
% move all data, so that d1 is in the origin (0,0)
dataNorm = bsxfun(@minus,D,D(d1,:));
%expand points by maxD
dataNorm = bsxfun(@times, dataNorm, maxD);
% move back to original position of d1
D = bsxfun(@plus,dataNorm,D(d1,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






