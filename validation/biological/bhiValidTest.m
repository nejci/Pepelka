
% export GO base to R list, delete trailing commas in output file before use.
% In R command window run: 
% source("D:\\Nejc\\LASPP\\mojeDelo\\KI\\Horvat\\matlab\\bhiValidTest_all.txt")

goNames=unique({GOstruct.functionID});
n = length(goNames);

fd = fopen(['bhiValidTest_',num2str(n),'.txt'],'w');
fprintf(fd,'goNames = c(');
fprintf(fd,'"%s", ',goNames{1:n});
fprintf(fd,')\n');

fprintf(fd,'genes = list(');
sumGenes = 0;
for i=1:n
    mask=strcmp(goNames{i},{GOstruct.functionID});
    sumGenes = sumGenes + sum(mask);
    c = {GOstruct(mask).geneID};
    fprintf(fd,'c("%s"', c{1});
    if length(c)>1
        fprintf(fd,', "%s"',c{2:end});
    end
    fprintf(fd,'),\n');
end  
fprintf(fd,')\n');
fprintf(fd,'names(genes) = goNames\n');

fclose(fd);
disp(sumGenes);