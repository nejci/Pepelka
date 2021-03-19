function BSI = indexBSI(genenames, labels, labelsDel, annotations, aspect, ignoreEvidence)

% Function computes Biological Stability Index (Datta and Datta 2006)
%-------------------------------------------------------------------------
% INPUTS
%   genenames       (cell)     strings with gene identifiers (must match 
%                              with identifiers in annotations)
%
%   labels          (vector)   cluster labels
%
%	labelsDel		(matrix)	matrix [n X d] of labels of data when 
%								clustering data without corresponding
%								column, i.e. labelsDel(:,i) contains labels
%								when clustering data without column i.
%
%   annotations     (struct)   structure array with fields: 
%                               .geneID         (gene identifier, e.g., Affymetrix probe ID)
%                               .functionID     (functional class identifier, e.g., GO ID)
%                               .aspect         (optional - when used with GO)
%                               .evidence       (optional - when used with GO)
%
%   aspect          (cell)      which aspects to include:
%                               'BP' | 'CC' | 'MF' (or any subset of them); 
%                               leave empty to compute all of them.
%
%   ignoreEvidence  (cell)      which evidence to ignore evidence code,
%                               i.e.: 'EXP', 'IDA', 'IEA',... or any subset;
%                               leave empty to compute all of them.
%--------------------------------------------------------------------------
% OUTPUTS:
%   value			(scalar)	value of BSI index
%
%--------------------------------------------------------------------------
% REQUIRES:                     Statistics Toolbox (crosstab)
%
%------- LEGAL NOTICE -----------------------------------------------------
% Copyright (C) 2012  Nejc Ilc
% Part of Pepelka package based on R package clValid.
%
%------- REFERENCE --------------------------------------------------------
% Datta, S., & Datta, S. (2006). 
% Methods for evaluating clustering algorithms for gene expression data 
% using a reference set of functional classes. 
% BMC bioinformatics, 7, 397. doi:10.1186/1471-2105-7-397
%
%------- VERSION ----------------------------------------------------------
% Version: 1.0
% Last modified: 6-Nov-2012 by Nejc Ilc
%
%------- CONTACT ----------------------------------------------------------
% Please write to: Nejc Ilc <myName.mySurname@fri.uni-lj.si>
%==========================================================================

% tic
if exist('ignoreEvidence','var') && ~isempty(ignoreEvidence)
    evidenceMask = ~ismember({annotations.evidence}, ignoreEvidence);
else
    evidenceMask = 1;
end

if exist('aspect','var') && ~isempty(aspect)
    aspectMask = ismember({annotations.aspect}, aspect);
else
    aspectMask = 1;
end

geneMask = ismember({annotations.geneID}, genenames);

mask = (geneMask & evidenceMask & aspectMask);
geneID = {annotations(mask).geneID};
functionID = {annotations(mask).functionID};


[FC,ia,ic] = unique(functionID); % function classes
nFC = length(FC);

% prepare data once to speedup for-loops
% For every function class, gather genes (their indeces) that belong to it

goList = cell(nFC,1);
nAnnot = zeros(nFC,1);
for g=1:nFC
    FCgenes = unique(geneID(ic == g)) ;
    nAnnot(g)=length(FCgenes);
    goList{g} = find(ismember(genenames, FCgenes));
end
% toc


BSIvec = zeros(1,nFC);

[n,d] = size(labelsDel);

BSI = zeros(1,d);

% loop over deleted columns
for del = 1:d
    labelsDel_i = labelsDel(:,del);
    
    overlap = crosstab(labels,labelsDel_i);
    rsums = sum(overlap,2);
    
    for FCk = 1:nFC
       osum = 0;
       gIdx = goList{FCk}; %find which genes are annotated
       if (length(gIdx)<2) % only one gene, skip
           continue;
       end
       %all pairs for labels->labelsDel
       n_gIdx = length(gIdx);
       for gxI = 1:(n_gIdx-1)
          for gyI = (gxI+1) : n_gIdx
            
              i = labels(gIdx(gxI));
              j = labelsDel_i(gIdx(gyI));
              osum = osum + overlap(i,j)/rsums(i);
            
          end
       end
       %all pairs for labelsDel->labels
       for gxI = 1:(n_gIdx-1)
          for gyI = (gxI+1) : n_gIdx
            
              i = labels(gIdx(gyI));
              j = labelsDel_i(gIdx(gxI));
              osum = osum + overlap(i,j)/rsums(i);
            
          end
       end
       BSIvec(FCk) = osum/(nAnnot(FCk)*( max(nAnnot(FCk)-1,1) ));
    end
    BSI(del) = mean(BSIvec); 
    
end

BSI = mean(BSI);
% toc
