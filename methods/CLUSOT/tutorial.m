function tutorial
% 2d map for the Iris data set
cd(clusot_root);
load iris.mat;
sM=som_make(iris_train);
plot_top(sM,iris_train);
[clust,sf,sl,C,BB,Mx,pMx,BP]=clusot2d(sM,iris_train);
plot_clusot(sf,sl,C,BP,pMx);
[clust,sf,sl,C,BB,Mx,pMx,BP]=clusot2d(sM,iris_train,'theta',0.5);
plot_clusot(sf,sl,C,BP,pMx);
plot_hit_num(sM,iris_train);
sM=som_make(iris_train,'lattice','rect');
plot_top(sM,iris_train);
[clust,sf,sl,C,BB,Mx,pMx,BP]=clusot2d(sM,iris_train);
plot_clusot(sf,sl,C,BP,pMx);

% 3D map for the Wine data set
cd(clusot_root);
load wine.mat
sM=som_make(wine_train,'lattice','rect','msize',[3 4 5]);
[clust,sf,sl,C,Mx,pMx,BP]=clusotnd(sM,wine_train);
plot_clusot3d(sM,wine_train,'C',C,'clust',clust,'pMx',pMx);

% Using external SOM data
cd(clusot_root);
[sM,f]=kd2som('11x9_map.txt', '11x9_hits.txt');
plot_top(sM,f)
[clust,sf,sl,C,BB,Mx,pMx,BP]=clusot2d(sM,f);
plot_clusot(sf,sl,C,BP,pMx);
[clust,sf,sl,C,Mx,pMx,BP]=clusotnd(sM,f,'theta',0.5);
plot_clusot(sf,sl,C,BP,pMx)

end % tutorial
