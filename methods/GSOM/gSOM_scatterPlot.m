function [fig]=gSOM_scatterPlot(data,labels,K,options)

% [fig]=plot_scatterPlot(data,labels,K,options)
%--------------------------------------------------------------------------
% Function pplk_scatterPlot plots a scatter plot of
% data clustering performed on data resultng in vector labels. 
%--------------------------------------------------------------------------
% INPUTS		
%	data	: (matrix) input data, size of [N x D]
%	labels	: (vector) clustering result
%	K		: (int)	number of clusters
%	options	: (struct)	options structure; can be
%						[] or non-existent for defaults;
%						Fields:
%						.title
%						.axisLabels (1x2 cell of strings)
%						.axisTicks (1x2 cell of vectors)
%						.colorMode ['color'|'pattern'|'mixed']
%						.normalize [0|1] - normalize on interval [0,1]
%						.fig - force this figure handle
%--------------------------------------------------------------------------
% OUTPUTS		
%	fig		: (handle)	handle on figure
%--------------------------------------------------------------------------
% Copyright (C) 2011  Nejc Ilc
%==========================================================================

%DEFAULTS
axisTitle='';
axisLabels={'',''};
axisTicks={[],[]};
colorMode='color';
normalize=0;
		
if(exist('options','var'))
	if isfield(options,'title')
		axisTitle=options.title;
	end
	if isfield(options,'axisLabels')
		axisLabels=options.axisLabels;
	end
	if isfield(options,'axisTicks')
		axisTicks=options.axisTicks;
	end
	if isfield(options,'colorMode')
		colorMode=options.colorMode;
	end
	if isfield(options,'normalize')
		normalize=options.normalize;
	end
	if isfield(options,'fig')
		figure(options.fig);
	else
		fig=figure();
	end
end

[N,D]=size(data);

if normalize
    data=pplk_normalize(data,'range');
end

%marker size
mS=20;

if(K<5)
    colors=[0 0 153/255; 0 153/255 1; 153/255 1 102/255; 1 102/255 0];

else
    colors=colormap(jet(K*10));
    colors=colors(1:10:end,:);
end

if(D>2)
    %PCA on 2D
    n=2;
    eig_vec=princomp(data,'econ');
    
    newData=(eig_vec(:,1:n)' * data');
    if normalize
        data=pplk_normalize(newData','range');
    else
        data=newData';
    end
end

patterns={'*','x','o','square','diamond','+','.','^','v','>','<','pentagram','hexagram'};
nP=numel(patterns);


%COLOR MODE - color
if(strcmp(colorMode,'color'))
    
	for i=1:K
       plot(data(labels==i,1),data(labels==i,2),'.','color',colors(i,:),'markersize',mS);
       hold on;
	end
    hold off;
    

%COLOR MODE - patterns
elseif (strcmp(colorMode,'pattern'))
    %marker size
	mS=10;
	
	for i=1:K
       plot(data(labels==i,1),data(labels==i,2),['k',patterns{mod(i-1,nP)+1}],'markersize',mS);
       hold on;
	end
	hold off;
	
   
%COLOR MODE - mixed (colored patterns)
else
	%marker size
	mS=10;
	
	for i=1:K
       plot(data(labels==i,1),data(labels==i,2),patterns{mod(i-1,nP)+1},'color',colors(i,:),'markersize',mS);
       hold on;
	end
	hold off;
end

if ~isempty(axisTicks{1}) && ~isempty(axisTicks{2})
	set(gca, 'XTick',axisTicks{1},'YTick',axisTicks{2});
end

title(axisTitle);
xlabel(axisLabels{1});
ylabel(axisLabels{2});

axis('equal');
axis('square');

end