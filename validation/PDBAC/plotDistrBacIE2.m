%This function plot the posterior distribution of the balanced accuraccy of a
%confusion matrix from a binary or multi-class classifer.
%In the figure it is also indicated the chance's treshold, mean, mode, median an
%alpha-confidence interval.
%
% Usage:
%     plotDistrBacIE2(AB,alpha,bmean,CI,chanceP,bmode,bdist,res,loc_leg)
%
% Modified by: Nejc Ilc
%
% Edited by:
% Henry Carrillo, University of Zaragoza, Spain
% http://www.hcarrillo.com/
%
% Original coder:
% Kay H. Brodersen, ETH Zurich, Switzerland
% http://people.inf.ethz.ch/bkay/
% -------------------------------------------------------------------------
function plotDistrBacIE2(AB,alpha,bmean,CI,chanceP,bmode,bdist,res,loc_leg)

b_lower = CI(1);
b_upper = CI(2);

x = 0:res:1;
% Plot balanced accuracy
hold on;
inner = (b_lower <= x) & (x <= b_upper);
plotfill(x(inner), bdist(inner), [0.6 0.6 0.6]);
plot([chanceP chanceP], [0 betaavgpdf(bmode, AB, res)], '--', 'color', [0 0 0], 'linewidth', 2);
plot([bmean bmean], [0 betaavgpdf(bmean, AB, res)], 'color', [192,0,0]/255, 'linewidth', 2);
plot(x,bdist,'k');
xlabel('$$\lambda $$','interpreter','latex');
ylabel('$$p  \left( \lambda \vert \mathcal{D} \right) $$','interpreter','latex')
hleg1 = legend([num2str(round((1-alpha)*100)), '% CI = ', '[',num2str(b_lower,'%6.2f'),', ',num2str(b_upper,'%6.2f'),']'], ['chance = ', num2str(chanceP*100,'%6.1f'),'%'],...
    ['mean = ', num2str(bmean*100,'%6.1f'),'%'], 'Location', loc_leg);%'NorthWest' SouthEast SouthWest

set(hleg1,'FontSize',12);
v = axis;
v(3)=0.0;
axis(v);
hold off;
end