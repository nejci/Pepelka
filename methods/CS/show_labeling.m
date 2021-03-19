% Plot data set  
function show_labeling(data,labels,f);
   if (size(data,1) > 2);
      data = pca(data,2);
   end;
   directed_graph_labels(data,labels,f)
   pause(1)
