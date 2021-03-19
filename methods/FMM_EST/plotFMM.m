function h = plotFMM(data,probs,means,covars)

parentDir = chdir('online');

PlotData(data');
PlotGM(probs,means,covars);

chdir(parentDir);