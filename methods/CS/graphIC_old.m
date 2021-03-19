%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB function: graphIC.m 
% Author: Robert Jenssen
% Last update: November 6, 2008. Nejc Ilc.
%               
%               - dodano risanje po korakih
%               - dodano odkrivanje kolena funkcije IC - pravega stevila
%               gruè
% 
% NOTE: Clustering algorithm implementing the 
% method described in R. Jenssen, J. C. Principe and 
% T. Eltoft, "Cauchy-Schwartz pdf Divergence Measure for 
% non-Parametric Clustering," in  NORSIG2003, Bergen, 
% Norway, October 2003. Recommended to use sigh=sigl,
% since it was just an experiment trying to see whether an 
% annealing of the kernel size would make sense. 
%
% OUTPUT
% l:    Vector containing the label for each data pattern (for K_true)
% IC:   Vector containing value of IC for every K (Kin:Kf)
% l_h:  Matrix containing history of labels
%
% INPUT
% data: Data set to be clustered
% sigh: Start kernel size - sigma high
% sigl: Stop kernel size - sigma low
% rate: Kernel size decrease
% Kf:   Final number of clusters
% Kin:  Initial number of clusters
% Nin:  Initial size of the Kin clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [l,K_true,IC,l_h]=graphIC_old(data,sigh,sigl,rate,Kf,Kin,Nin,draw);

[dim,n]=size(data);

%history of labels
l_h=zeros(Kin-Kf+1,n);

if(draw & dim==2)
   fig=figure();
   fig2=figure();
   axis=[-1 1 -1 1];
   figure(fig);
   plot(data(1,:),data(2,:),'o');
   
   pause;
   
   % Seed the Kin clusters; return l (labels) and M (means) 
    [l,ul,KM]=seed(data,sigh,Kin,Nin,fig2);
    
else
    % Seed the Kin clusters; return l (labels) and M (means) 
[l,ul,KM]=seed(data,sigh,Kin,Nin,0);
end



if(draw & dim==2)
   directed_graph_labels(data,l,fig2);
   pause;
end

% Determine volume and cut contributions of each cluster
[VOL,CUT]=cut_vol(l,data,sigh,Kin);


% Current number of clusters is Kin and unlabeled points
K = Kin; unlabeled = length(find(ul == 1));

it = 1;
% Loop until final number of clusters is reached
while (K >= Kf);
   
    % Exponential decrease in kernel size
    sig = sigl + (sigh-sigl)*exp(-rate*(it-1));
    
    % Inner loop until no more unlabeled points
    while (unlabeled > 0);
        
        % Find next xi to be clustered (nearest some mean)
        %d[xi_index]=next_xi_mean(KM,data,ul,sig);
        
        % Find next xi to be clustered (nearest some unlabeled)
        [xi_index]=next_xi_ul(data,ul,sig);
        
        %plotit
        
        % Cluster xi according to IC
        [l,ul,KM,CUT,VOL]=clusterIC(l,ul,xi_index,data,KM,CUT,VOL,sig,K);
       
        if (draw==2 & dim==2)
        % Show adding of every pattern
        
            directed_graph_labels(data,l,fig2);
            title(['Dodajanje vzorca: ', int2str(xi_index) ]);
            pause
        end
        
        % One less unlabeled point
        unlabeled = unlabeled - 1;
        
    end; % while (unlabeled == 0);
    
    if (draw==1 & dim==2)
        
        
            directed_graph_labels(data,l,fig2);
            title('Dodajanje neoznacenih');
            pause
    end
    
    % Record value of IC
    %fCUT(it) = (0.5*sum(sum(CUT))); fVOL(it) = sqrt(prod(VOL));
    IC(it) = (0.5*sum(sum(CUT))) / sqrt(prod(VOL));
    l_h(it,:)=l;
    
    %if (K==4);
    %    l4 = l;
    %end;
    
    % Only if further re-labeling is to be done
    if (K > Kf);
        
        % Evaluate IC; set labels of 'worst' cluster to 0
        [l,ul,KM,CUT,VOL]=evaluateIC(l,ul,CUT,VOL,KM,K);
    
        if (draw & dim==2)
        % Show it!
        
            directed_graph_labels(data,l,fig2);
            title('Najslabsa skupina');
            pause
        end
        
        % Need to update unlabeled;
        unlabeled = length(find(ul == 1));
        
    end; % if (K > Kf);
      
    % Number of clusters reduced by one
    K = K - 1;
    
    % Iterate
    it = it + 1;
        
end; % while (K > Kf);

%Added: we want to determine true number of clusters - where IC is minimal
%[IC_val,IC_ind]=min(IC);

% %prvi odvod
% diff=IC(2:end)-IC(1:end-1);
% %drugi odvod
% diff2=diff(2:end)-diff(1:end-1);
% 
% [val,ind]=max(diff2);
% 
% K_true=Kin-ind;
% l=l_h(ind+1,:);
% 
% if(draw & dim==2)
%    
%     directed_graph_labels(data,l,fig2);
%     title(['Konec, K=',int2str(K_true)]);
%     
%     figure();
%     plot(IC);
%     hold on;
%     
%     %figure();
%     plot(1.5:1:(Kin-Kf+1-0.5),diff,'r.');
%     
%     %figure();
%     plot(2:1:(Kin-Kf),diff2,'g*');
%     
% end

l=l_h(end,:);
K_true=K;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions used by graphIC.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Seeds the initial clusters and calculates mean 
function [l,ul,KM]=seed(data,sig,Kin,Nin,fig);

    [d,n]=size(data);

    % No points are labeled initially (0)
    l = zeros(1,n);
    
    %Dal sem ga izven zanke, da se lahko kontrolira
    R = randperm(n);
    %R=[1,4,6,8]; %za test
    
    for i = 1 : Kin;
        
        % Create a bit string 'ul' with 1 = unlabeled
        ul = zeros(1,n); f0 = find(l == 0); 
        ul(f0) = ones(1,length(f0));
        
        % Random unlabeled seed point si
        stop = 0; it = 1;
        while (stop == 0);
            si = R(it);
            if (ul(si) == 1);
                stop = 1;
            else;
                it = it + 1;
            end;
        end;
       
        % si'th data vector
        xsi = data(:,si);
        
        if(fig>0)
           figure(fig);
            plot(xsi(1), xsi(2), 'r+', 'markersize', 10);
            hold on;
        end
        
        % Similarity between si and xj's
        [Gsi_xj]=Gijs(data,xsi,sig);
        
        % Set to zero distance to labeled points; sort
        [sorted,index]=sort(Gsi_xj .* ul);
        
        % Sort in descending order
        index = fliplr(index);      
          
        % Label Nin points nearest to si with label i - oznacevanje Nin
        % najblizjih zacetnemu - morebitna izboljšava: najbljižji
        % kateremukoli v skupini
        l(index(1:Nin)) = i * ones(1,Nin); 
        
        % Update unlabeled vector ul
        ul = zeros(1,n); f0 = find(l == 0); 
        ul(f0) = ones(1,length(f0));
        
        % Calculate mean of cluster i
        KM(:,i) = mean(data(:,index(1:Nin)),2); % mean
        
    end; % for i = 1 : Kin; 
    
% Determine volume and cut contributions of each initial cluster
function [VOL,CUT]=cut_vol(l,data,sig,Kin);
    
    % Initialize
    VOL = zeros(1,Kin); CUT = zeros(Kin);

    % Loop over all Kin clusters
    for i = 1 : Kin;
        
        % Find index of members of cluster i
        fi = find(l==i);
        
        % Second loop over all Kin clusters
        for j = 1 : Kin;
            
            % Members of cluster j
            fj = find(l==j);
            
            % CUT contribution; cluster i --> cluster j - števec v formuli
            % (str. 2)
            if (i~=j);
                
                % Do the same for all members of cluster i
                for k = 1 : length(fi);
                    
                    % Member fi(k) of cluster i
                    xk = data(:,fi(k));
                    
                    % Find similarity between xk and members of j 
                    [Gxk_dataj]=Gijs(data(:,fj),xk,sig);
                
                    % Update CUT(i,j)
                    CUT(i,j) = CUT(i,j) + sum(Gxk_dataj);
                    
                end; % for k = 1 : length(fi);
            
            % Volume of cluster i - imenovalec v formuli
            else;
                
                % Do the same for all members of cluster i
                for k = 1 : length(fi);
                    
                    % Member fi(k) of cluster i
                    xk = data(:,fi(k));
                    
                    % Find similarity between xk and members of j=i 
                    [Gxk_dataj]=Gijs(data(:,fj),xk,sig);
                
                    % Update VOL(i)
                    VOL(i) = VOL(i) + sum(Gxk_dataj);
                    
                end; % for k = 1 : length(fi);
                
            end; % if (i~=j);
        end; % for j = 1 : Kin;
    end; % for i = 1 : Kin;
                
% Find next xi to be clustered (nearest some mean)
function [xi_index]=next_xi_mean(KM,data,ul,sig);
    
    % Current number of clusters
    K = size(KM,2);
    
    % Loop over all cluster means
    for i = 1 : K;
        
        % i'th Mean
        Mi = KM(:,i);
        
        % Similarity; Mi to all xj in data set
        [GMi_xj]=Gijs(data,Mi,sig);
        
        % Sort to find the closest (ascending order)
        [sorted,index]=sort(GMi_xj .* ul);
        
        % Store largest similarity value and index
        NN(1,i) = sorted(end);
        NN(2,i) = index(end);
        
    end; % for i = 1 : K;
    
    % Find index xi overall most similar to some Mi
    [xi_sorted,xi_index]=sort(NN(1,:));
    xi_index = NN(2,xi_index(end));
    
% Find next xi to be clustered (nearest to some already labeled)
function [xi_index]=next_xi_ul(data,ul,sig);
   
    % Find the number of unlabeled data patterns and members
    ful = find(ul == 1);
    dataul = data(:,ful);

    % Loop over all unlabeled points
    for i = 1 : length(ful);
        
        % i'th unlabeled point
        uli = dataul(:,i);
        
        % Similarity; uli to all xj in data set
        [Guli_xj]=Gijs(data,uli,sig);
        
        % Sort to find the closest to some labeled 
        [sorted,index]=sort(Guli_xj .* (-1*(ul-1)));
        
        % Store largest similarity value and index
        NN(1,i) = sorted(end);
        NN(2,i) = index(end);
        
    end; % for i = 1 : K;
    
    [maxGij,i_index]=max(NN(1,:));
    i_index = i_index(1);
    xi_index = ful(i_index);
    
    
% Cluster xi according to IC
function [l,ul,KM,CUT,VOL]=clusterIC(l,ul,xi_index,data,KM,CUT,VOL,sig,K);

    % Data pattern pointed to by xi_index
    xi = data(:,xi_index);

    % Similarity of xi to itself - popravljeno - dodana 2 
    [d,n] = size(data);
    Gxi_xi = inv((2*pi)^(d/2)*sig^d);
    
    % Overall similarity from xi to each of K clusters
    for j = 1 : K;
        
        % Index of members of cluster j
        fj = find(l == j);
        
        % Members of cluster j
        dataj = data(:,fj);
        
        % Kernel similarity from xi to members of cluster j
        [Gxi_dataj]=Gijs(dataj,xi,sig);
         
        % Overall similarity from xi to members of cluster j
        SimGxi_dataj(j) = sum(Gxi_dataj);
        
    end; % for j = 1 : K;
    
    % Calculate IC for xi assigned to each of K clusters
    for j = 1 : K;
        
        % Assign xi to cluster j, adjust volume of cluster j  
        tmpVOL = VOL;
        tmpVOL(j) = VOL(j) + SimGxi_dataj(j) + Gxi_xi;
    
        % tmpSimGxi_dataj used to calculate CUT
        tmpSimGxi_dataj = SimGxi_dataj;
        tmpSimGxi_dataj(j) = 0; %ker so x_i in ostali iz clustra j v istem clustru: M_ij=0
        
        % Temporary matrix of CUT contributions
        tmpCUT = CUT;
        tmpCUT(:,j) = tmpCUT(:,j) + tmpSimGxi_dataj';
        tmpCUT(j,:) = tmpCUT(j,:) + tmpSimGxi_dataj;
        
        % CUT value with xi assigned to cluster j
        CUTxi_to_j = sum(sum(tmpCUT));
        
        % IC with xi assigned to cluster j
        IC(j) = (0.5 * CUTxi_to_j) / sqrt(prod(tmpVOL)); 
        
    end; % for j = 1 : K;
    
    % Determine which cluster xi should be assigned to
    [minIC,j] = min(IC); minIC = minIC(1); j = j(1);
    
    % New volume of cluster j with xi counted!
    VOL(j) = VOL(j) + SimGxi_dataj(j) + Gxi_xi;
    
    % Update CUT contributions!
    tmpSimGxi_dataj = SimGxi_dataj;
    tmpSimGxi_dataj(j) = 0;
    CUT(:,j) = CUT(:,j) + tmpSimGxi_dataj';
    CUT(j,:) = CUT(j,:) + tmpSimGxi_dataj;
    
    % Update label
    l(xi_index) = j;
    
    % Update vector over unlabeled patterns
    ul(xi_index) = 0;
    
    % Brute force update of new mean
    fj = find(l == j);
    KM(:,j) = mean(data(:,fj),2); 
        
% Kernel evaluated distance from xi to all other xj    
function [Gij]=Gijs(data,xi,sig);
    
    % Number and dimension of data vectors 
    [d,n] = size(data);
    %mislim, da manjka 2 pri sig^2, ker pride iz tega sqrt(det(2*sigma^2*I))
    C1 = inv((2*pi)^(d/2)*sig^d);
    %mislim, da manjka še ena 2 pri 2*sig^2
    C2 = -inv(2*sig^2);

    Xi = zeros(d,n);
    for t = 1:d
        Xi(t,:) = xi(t)*ones(1,n);
    end;
    
    diff =  Xi - data;
    
    % Kernel evaluations
    Gij = C1*exp(C2*sum(diff.*diff,1));     
     
% Evaluate IC; set labels of 'worst' cluster to 0
function [l,ul,KM,CUT,VOL,minIC]=evaluateIC(l,ul,CUT,VOL,KM,K);

    % Calculate IC resulting from removing cluster i
    for i = 1 : K;
        
        % Remove volume of cluster i
        tmpVOL = VOL;
        tmpVOL(i) = 1;

        % Remove CUT contribution from cluster i
        tmpCUT = CUT;
        tmpCUT(:,i) = [];
        tmpCUT(i,:) = [];
        
        % Calculate IC with cluster i removed
        IC(i) = (0.5*sum(sum(tmpCUT)))/sqrt(prod(tmpVOL));

    end; % for i = 1 : K;
    
    % Identify cluster which makes IC smallest (largest CS divergence) when removed
    [minIC,index_minIC]=min(IC); 
    minIC=minIC(1); 
    index_minIC=index_minIC(1);
    
    % Update VOL
    VOL(index_minIC) = [];
    
    % Update CUT
    CUT(:,index_minIC) = [];
    CUT(index_minIC,:) = [];
    
    % Set labels of cluster i to 0 (free to be re-clustered)
    findex_minIC = find(l==index_minIC);
    l(findex_minIC) = zeros(1,length(findex_minIC));
    
    % Update the rest of the labels (to a value one less than before)
    f_adjust = find(l > index_minIC);
    l(f_adjust) = l(f_adjust) - 1;
    
    % Create a bit string 'ul' with 1 = unlabeled
    ul = zeros(1,length(l)); f0 = find(l == 0); 
    ul(f0) = ones(1,length(f0));
    
    % Remove mean corresponding to what was cluster index_minIC
    KM(:,index_minIC) = [];
