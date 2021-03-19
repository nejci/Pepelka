Code for Bayesian Cluster Ensembles


The folder is the matlab code for Bayesian Cluster Ensembles (BCE). The model was proposed in the paper: H. Wang, H. Shan, A. Banerjee. Bayesian Cluster Ensembles. Statistical Analysis and Data Mining, 2011.

The zip file contains the following files:

runBCE.m:               An example on how to run the code.
learnBCE.m:             Learn BCE. It calls bceEstep.m and bceMstep.m.
bceEstep.m:             Variational E-step.
bceMstep.m:             Variational M-step.
calculateAccuracy.m     Compute the clustering accuracy.
Iris.mat:               Sample data.
readme.txt:             Readme file.