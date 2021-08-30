function [fig,data] = pplk_scatterPlot(data,labels,K,options)
% [fig,data]=PPLK_SCATTERPLOT(data,labels,K,options)
% Function plots a scatter plot of data and clustering results performed on
% data resultng in a labels vector.
%
% INPUTS		
%   data
%       A N-by-D matrix of input data.
%
%   labels
%       - A vector of data labels - clustering.
%       - A matrix set of clusterings - plot them on the subfigure.
%
%   K		
%       Number of clusters; if [] or non-existent, K becomes the number of 
%       unique labels.
%
%   options	
%       Options structure; can be [] or non-existent for defaults.
%       Fields:
%           .title 
%               Main title.
%           .subtitle 
%               Cell of strings - one string for each subplot.
%           .axisLabels 
%               1-by-2 cell of strings.
%           .axisTicks 
%               1-by-2 cell of vectors.
%           .axisStyle 
%               String, eg. 'square', 'equal', 'tight'.
%           .colorMode 
%               ['color'|'pattern'|'mixed'].
%           .normalize 
%               Normalize on interval [0,1].
%           .fig 
%               Force this figure handle.
%           .annotations 
%               n-by-1 cell of strings, where n is number of data points.
%           .markerSize
%               Size of the marker.
%           .interpreterMode 
%               'tex','latex',['none'].
%
%
% OUTPUTS		
%   fig	
%       Handle on figure.
%
%   data
%       Matrix of processed data displayed on graph (normalized, PCA 
%       transformed if required).
%
%
% This is a part of the Pepelka package.
% Contact: Nejc Ilc (nejc.ilc@fri.uni-lj.si)
% https://github.com/nejci/Pepelka

[N,D]=size(data);

%DEFAULTS
axisTitle='';
axisSubtitle='';
axisLabels={'',''};
axisTicks={[],[]};
axisStyle = 'normal';
colorMode='color';
annotations={};
normalize = 0;
mS=20; % default marker size
interpreterMode = 'none'; %tex

if ~exist('K','var') || isempty(K)
    hasK = 0;
else
    hasK = 1;
end

if ~exist('labels','var') || isempty(labels)
    labels = ones(N,1);
end


if(exist('options','var'))
    if isfield(options,'title')
        axisTitle=options.title;
    end
    if isfield(options,'subtitle')
        axisSubtitle=options.subtitle;
    end
    if isfield(options,'axisLabels')
        axisLabels=options.axisLabels;
    end
    if isfield(options,'axisTicks')
        axisTicks=options.axisTicks;
    end
    if isfield(options,'axisStyle')
        axisStyle=options.axisStyle;
    end
    if isfield(options,'colorMode')
        colorMode=options.colorMode;
    end
    if isfield(options,'normalize')
        normalize=options.normalize;
    end
    if isfield(options,'fig')
        %figure(options.fig);
        fig=options.fig;
    else
        fig=figure();
    end
    if isfield(options,'annotations')
        annotations=options.annotations;
    end
    if isfield(options,'markerSize')
        mS=options.markerSize;
    end
    if isfield(options,'interpreterMode')
        interpreterMode=options.interpreterMode;
    end
    
else
    fig=figure();
end

% force labels vector to be a column vector
if isvector(labels)
    labels = labels(:);
end

% Super title
%suptitle(axisTitle); % requires Bioinformatisc toolbox
annotation('textbox', [0 0.9 1 0.1], ...
    'String', axisTitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center',...
    'Interpreter','none',...
    'FontSize',14);

%-------------------------------------------------------------------------

switch normalize
    case {0,'none'}
        ;
    case {1,'range'}
        data=pplk_normalize(data,'range');
    case {2,'propor'}
        data=pplk_normalize(data,'propor');
    otherwise
        fprintf(1,'Wrong normalization method! No normalization will occur.\n');
end

if(D>2)
    %PCA projection onto 2D
    n=2;
    %[~, newData]=princomp(data,'econ');
    [~, data]=pca(data,'NumComponents',n);
    %data=newData(:,1:n);
    
end

patterns={'+','o','x','^','square','diamond','*','v','>','<','pentagram','hexagram','.'};
nP=numel(patterns);


% produce subplot for every clustering in columns of labels
ensembleSize = size(labels,2);

% max 4 subplots in a row of panel
nCols = min(ensembleSize,4);

nRows=floor(ensembleSize/nCols);
if rem(ensembleSize,nCols)~=0
    nRows=nRows+1;
end



labelsAll = labels;

for plotInd = 1:ensembleSize
    labels = labelsAll(:,plotInd);
    if ~hasK
        K = length(unique(labels));
    end
    
    if(K<5)
        colors=[    218/255 37/255 29/255; ...  % red
                    40/255 22/255 111/255; ...  % dark blue
                    132/255 194/255 37/255; ... % green
                    232/255 188/255 0/255];     % yellow
        
    else
        colors=colormap(jet(K*10));
        colors=colors(1:10:end,:);
    end
    
    if ensembleSize > 1
        subplot(nRows,nCols,plotInd);
    end
    %COLOR MODE - color
    if(strcmp(colorMode,'color'))
        
        for i=1:K
            plot(data(labels==i,1),data(labels==i,2),'.','color',colors(i,:),'markersize',mS);
            hold on;
        end
        hold off;
        
        
    %COLOR MODE - patterns
    elseif (strcmp(colorMode,'pattern'))
         
        for i=1:K
            plot(data(labels==i,1),data(labels==i,2),['k',patterns{mod(i-1,nP)+1}],'markersize',mS);
            hold on;
        end
        hold off;
        
        
    %COLOR MODE - mixed (filled, big patterns)
    elseif (strcmp(colorMode,'mixedClear'))
        patterns={'o','>','pentagram','v','square','diamond','^','hexagram','<','.','*','x','+'};
        
        
        for i=1:K
            plot(data(labels==i,1),data(labels==i,2),patterns{mod(i-1,nP)+1},'MarkerFaceColor',colors(i,:),'MarkerEdgeColor',colors(i,:),'markersize',mS);
            hold on;
        end
        hold off;
        
    else
        for i=1:K
            plot(data(labels==i,1),data(labels==i,2),patterns{mod(i-1,nP)+1},'color',colors(i,:),'markersize',mS);
            hold on;
        end
        hold off;
    end
    
    % plot annotations below data points as well
    if ~isempty(annotations)
        assert(length(annotations) == N);
        offset = (max(data(:,2))-min(data(:,2)))/100;
        for ind=1:N
            of = 0;%offset*2*rand();
            text(data(ind,1), data(ind,2)-of, annotations{ind},'FontSize',8,'Interpreter','None','HorizontalAlignment','center','VerticalAlignment','top','Rotation',0);
            
        end
    end
    
    if ~isempty(axisTicks{1}) && ~isempty(axisTicks{2})
        set(gca, 'XTick',axisTicks{1},'YTick',axisTicks{2});
    end
    
    if iscell(axisSubtitle)
        title(axisSubtitle(plotInd),'Interpreter',interpreterMode);
    end
    
    xlabel(axisLabels{1});
    ylabel(axisLabels{2});
    
    axis(axisStyle);
end

