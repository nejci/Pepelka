function [labelsCE]=runLCE(E, K, methodSim, methodCons, dc, iter)
%==========================================================================
% FUNCTION: labelsCE = runLCE(E, K, methodSim, methodCons, dc, iter)
% DESCRIPTION:	This function creates similarity matrix and performs 
%				the final clustering using Hierarchical algorithms 
%				(Single-Linkage:SL, Complete-Linkage:CL and
%				Average-Linkage:AL) as consensus methods.
%
% INPUTS:	E = an ensemble members matrix 
%			K = the prefered number of clusters
%			methodSim = name of the method for creating similarity matrix
%						('CTS','SRS','ASRS')
%			methodCons = name of the method for final clustering
%						 ('SL'-single, 'CL'-complete, 'AL'-average linkage)
%			dc = decay factor (scalar) - input parameter to methodSim [default=0.8]
%			iter = number of iterations for SRS method [default=5]
%
% OUTPUTS: labelsCE = an N-by-1 vector of clustering result 
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
% modified: Nejc Ilc
%==========================================================================


% default values for decay factor and num. of iterations (for SRS only)
if nargin<6; iter = 5; end;  
if nargin<5; dc = 0.8; end;   

methodSim=upper(methodSim);
methodCons=upper(methodCons);

if ~ismember(methodSim,{'CTS', 'SRS', 'ASRS'})
	error([methodSim, ' is NOT a valid method name! Use: ''CTS'', ''SRS'' or ''ASRS''.']);
end

if ~ismember(methodCons,{'SL', 'CL', 'AL'})
	error([methodCons, ' is NOT a valid method name! Use: ''SL'', ''CL'' or ''AL''.']);
end

if (dc > 1) || (dc < 0)
	error('LCE: "dc" must be in the range of 0 and 1.')
end


% Generating similarity matrices and performing consensus functions
%disp(['LCE: Generating ' upper(methodSim) ' matrix ...']);
switch methodSim
	case 'CTS'
		S = cts(E, dc);
	case 'ASRS'
		S = asrs(E, dc);
	case 'SRS'
		S = srs(E, dc, iter);
	otherwise
		error(['LCE: no similarity matrix generation function called: ',methodSim]);
end

% perform consensus functions - linkage methods
%disp(['LCE: Preforming ' upper(methodCons) ' clustering ...']);
labelsCE=consFuncLCE(S,K,methodCons);



function CR = consFuncLCE(S, K, method)
%==========================================================================
% FUNCTION: CR = consFuncLCE(S, K, method)
% DESCRIPTION: This function performs the final clustering using Hierarchical
%              algorithms (Single-Linkage:SL, Complete-Linkage:CL and
%              Average-Linkage:AL) as consensus methods.
%
% INPUTS:     S = an N-by-N similarity matrix
%             K = the prefered number of clusters
%
% OUTPUTS: CR = an N-by-1 matrix of clustering results from given method
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

d = stod(S); %convert similarity matrix to distance vector
switch method
	case 'SL'
		% single linkage
		Z = linkage(d,'single');
		CR = cluster(Z,'maxclust',K);
	
	case 'CL'
		% complete linkage
		Z = linkage(d,'complete');
		CR = cluster(Z,'maxclust',K);
	case 'AL'
		% average linkage
		Z = linkage(d,'average');
		CR = cluster(Z,'maxclust',K);
	otherwise
		error(['LCE: no final consensus function called: ',method]);
end