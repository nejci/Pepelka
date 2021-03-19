function directed_graph_labels(data,labels,f);

figure(f);
clf;

N = size(data,1);

numClusters=0;

for i=1:N
    if(sum(labels==i) > 0)
        numClusters=numClusters+1;
    end
end

colors=colormap(jet(numClusters*10));
colors=colors(1:10:end,:);

plot(data(:,1), data(:,2), 'ko', 'markersize', 5);
hold on;

for i=1:numClusters
   plot(data(labels==i,1),data(labels==i,2),'.','color',colors(i,:));
   hold on;
    
end

% for i = 1 : N;
%     if (labels(i)==1);
%         plot(data(1,i),data(2,i),'r*');
%         hold on;
%     end;
% end;
% 
% for i = 1 : N;
%     if (labels(i)==2);
%         plot(data(1,i),data(2,i),'bx');
%         hold on;
%     end;
% end;
% 
% for i = 1 : N;
%     if (labels(i)==3);
%         plot(data(1,i),data(2,i),'yo');
%         hold on;
%     end;
% end;
% 
% for i = 1 : N;
%     if (labels(i)==4);
%         plot(data(1,i),data(2,i),'gs');
%         hold on;
%     end;
% end;
% 
% for i = 1 : N;
%     if (labels(i)==5);
%         plot(data(1,i),data(2,i),'md');
%         hold on;
%     end;
% end;
% 
% for i = 1 : N;
%     if (labels(i)==6);
%         plot(data(1,i),data(2,i),'v');
%         hold on;
%     end;
% end;
% 
% for i = 1 : N;
%     if (labels(i)==7);
%         plot(data(1,i),data(2,i),'^');
%         hold on;
%     end;
% end;
% 
% for i = 1 : N;
%     if (labels(i)==8);
%         plot(data(1,i),data(2,i),'>');
%         hold on;
%     end;
% end;
% 
% for i = 1 : N;
%     if (labels(i)==9);
%         plot(data(1,i),data(2,i),'<');
%         hold on;
%     end;
% end;
% 
% for i = 1 : N;
%     if (labels(i)==0);
%         plot(data(1,i),data(2,i),'m+');
%         hold on;
%     end;
% end;