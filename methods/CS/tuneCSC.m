function [err_h, labels_h, interval, IC_h]=tuneCSC(data,mode,interval,repeat,other,target,K)
%data - dimensions x patterns matrix
%mode - 'sigma', 'Kin', 'Nin' choice
%interval - vector of values for tuning
%repeat - number of repeated tests for each value of observed variable
%other - vector of other values as follows:
%         'sigma' - other=[Kin, Nin]
%         'Kin'   - other=[sigma, Nin]   
%         'Nin'   - other=[sigma, Kin] 
%
%target - correct labels in a vector
%K - true number of clusters
%
%err_h - history of errors (r x steps matrix)


steps=length(interval);

err_h=zeros(repeat,steps);
IC_h=zeros(repeat,steps);
labels_h=zeros(repeat,steps,size(data,2));

switch mode
    
    case 'sigma'
       for r=1:repeat 
            for i=1:steps

                [l,K_true,IC,l_h]=graphIC_old(data,interval(i),interval(i),1,K,other(1),other(2),0);
                IC_h(r,i)=IC(end);
                err_h(r,i)=evaluate(l_h(end,:),target,K);
                labels_h(r,i,:)=l_h(end,:);
                fprintf('Napredek - r=%d, step=%d\n',r,i);
            end
            
       end
       
    case 'Kin'
        for r=1:repeat 
            for i=1:steps

                [l,K_true,IC,l_h]=graphIC_old(data,other(1),other(1),1,K,interval(i),other(2),0);
                err_h(r,i)=evaluate(l_h(end,:),target,K);
                IC_h(r,i)=IC(end);
                labels_h(r,i,:)=l_h(end,:);
                fprintf('Napredek - r=%d, step=%d\n',r,i);
            end
        end
       
    case 'Nin'
        for r=1:repeat 
            for i=1:steps

                [l,K_true,IC,l_h]=graphIC_old(data,other(1),other(1),1,K,other(2),interval(i),0);
                err_h(r,i)=evaluate(l_h(end,:),target,K);
                IC_h(r,i)=IC(end);
                labels_h(r,i,:)=l_h(end,:);
                fprintf('Napredek - r=%d, step=%d\n',r,i);
            end
       end
end

plot(interval,mean(err_h,1),'-xr');title('Mean error');
figure;
errorbar(interval,mean(err_h),std(err_h));title('Errorbar');