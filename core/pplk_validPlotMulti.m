function [fig] = pplk_validPlotMulti(validResults, numClust, validMethods, clustMethods)

% validResults  (cell)      cell of [numValidMethods X length(numClust)]
%                           matrices with values of indeces in rows
%
% numClust      (vector)    vector of integers; for each index value 
%                           (columns of validResults) it gives the number
%                           of clusters
%
% validMethods  (cell)      cell of validation methods used to validate
%                           clustering results
%
% clustMethods  (cell)      cell of clustering methods used to produce
%                           results

numClust_len=length(numClust);
numClustMethods = length(clustMethods);

if ~exist('validMethods','var') || isempty(validMethods)
        [numValidMethods, numClust_len]=size(validResults{1});        
        validMethods = cellstr([repmat('index ',numValidMethods,1),  num2str([1:numValidMethods]')]);
end

numValidMethods = length(validMethods);

indexMatrix = zeros(numClustMethods,numClust_len);
indexCell = cell(1,numValidMethods);

for vInd = 1:numValidMethods
    for cInd = 1 : numClustMethods
        indexMatrix(cInd,:) = validResults{cInd}(vInd,:);
    end
    indexCell{vInd} = indexMatrix;
end


% max 4 subplots in a row of panel
nCols = min(numValidMethods,4);

nRows=floor(numValidMethods/nCols);
if rem(numValidMethods,nCols)~=0
	nRows=nRows+1;
end

% 
figure();
set(0, 'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1]);
%set(0, 'DefaultAxesColorOrder',[0 0 0]);
set(0, 'DefaultAxesLineStyleOrder','-|--|:|-.');


for validMethodInd = 1:numValidMethods


    subplot(nRows,nCols,validMethodInd);       
    plot(numClust, indexCell{validMethodInd}');

%     ht=title(['Internal validation']);
%     set(ht,'Interpreter','none');
    hl=legend(clustMethods,'Location','NorthEast');	
    set(hl,'Interpreter','none');
   
    if(length(numClust) > 1)
        xlim([min(numClust),max(numClust)]);
    end
    xlabel('k');
    ylabel(validMethods{validMethodInd});
    
end