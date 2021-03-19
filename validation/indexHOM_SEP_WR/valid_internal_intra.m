function [Hom,Sep,Cindex,Wintp] = valid_internal_intra(Smatrix,labels,dtype,id)
% indices base on intra and inter similarity

U= ind2cluster(labels);

Hom = 0;
Sep = 0;
Wintp = 0;
Dunn = 0;
DB = 0;
k = length(U);
nrow = size(Smatrix,1);

[avintra, avinter, intra, inter, nintra, ninter] = valid_intrainter(Smatrix,U);

ns = nrow*(nrow-1);         % number of all pairs of elements
for i = 1:nrow
    Smatrix(i,i) = 0;
end
cut = (sum(sum(Smatrix))-nrow)/ns;
ns = sum(nintra);
S = sum(intra);
[wintra, dbs] = find_nearpoint(Smatrix, cut, ns, 2);
for i = 1:nrow
    Smatrix(i,i) = inf;
end
[sinter, ni] = find_nearpoint(Smatrix, cut, ns, 1);
if ni < dbs
    wintra = wintra(1:ni);
elseif ni > dbs
    sinter = sinter(1:dbs);
end
R = Smatrix(wintra);
Smax = sum(R);
R = Smatrix(sinter);
Smin = sum(R);
Cindex = (S-Smin)/(Smax-Smin);        % C index, Hubert and Levin

if id
    return; % return when id==10, i.e., when C index is computed
end

ns = k*(k-1)/2;                                    % mates of clusters
[avintra, avinter, intra, inter, nintra, ninter] = valid_intrainter(1-Smatrix,U);
% Homogeneity & Separation
Hom = sum(intra)/sum(nintra);          % average intra similarity
Sep = Hom;
if ns > 0
    Sep = (sum(sum(inter)))/(sum(sum(ninter)));    % average inter
end
if strcmp(dtype,'euclidean') || (isnumeric(dtype) && dtype == 1)
    Hom = 1-Hom;
    Sep = 1-Sep;
elseif strcmp(dtype,'correlation') || (isnumeric(dtype) && dtype == 2)
    Hom = Hom+Hom -1;
    Sep = Sep+Sep -1;
end

% weight_intertra
wintra = zeros(1,k);
Sinter = zeros(1,k);
sinter = zeros(1,k);
Inter = inter+inter';
for i = 1:k
    ind = U{i};
    ni = length(ind);
    Sinter(i) = sum(Inter(i,:))/(nrow-ni);
    sinter(i) = sum(inter(i,:))/(nrow-ni);
    if ni ==1
        ni = 2;
    end
    wintra(i)=2*intra(i)/(ni-1);
end
if k == 2
    sinter(2) = 0.5*Sinter(2);
end
Sintra = sum(wintra);
%Sinter = sum(Sinter);
Sinter = sum(sinter);
Wint = 1-Sinter/Sintra;
Wintp = (1-2*k/nrow)*Wint;    % penalized

end

function [avintra, avinter, intra, inter, nintra, ninter] = valid_intrainter(Smatrix,U)
%caculate intra similarity/distance and inter similarity

NC=length(U);
nintra=zeros(1,NC);
intra=zeros(1,NC);
avintra=zeros(1,NC);
ninter=zeros(NC,NC);
avinter=zeros(NC,NC);

inter=zeros(NC,NC);
for i=1:NC
    ind=U{i};
    ni=length(ind);
    R=Smatrix(ind,ind);
    if ni==1
        nintra(i)=1;
        intra(i)=0.2;                 %single element similarity:?max(max(Smatrix))/2
        avintra(i)=0.2;
    else
        Q=sum(sum(triu(R,1)));
        T=sum(sum(tril(R,-1)));
        if Q < T
            Q = T;
        end
        intra(i)=Q;                             %similarity sum of intra cluster i
        if ni < 2
            nintra(i) = 1;
        else
            nintra(i)=(ni*(ni-1))/2;        %number of mates in half matrix
        end
        avintra(i)=intra(i)/nintra(i);     %average similarity in cluster i
    end
    
    for j=i+1:NC
        indj=U{j};
        nj=length(indj);
        R=Smatrix(ind,indj);
        inter(i,j)=sum(sum(R));          %disimilarity between cluster i & j
        ninter(i,j)=ni*nj;                       %number of non mates between
        if ninter(i,j) == 0
            ninter(i,j) = 1;
        end
        avinter(i,j)=inter(i,j)/ninter(i,j); %average similarity
    end
    
end

end

function [cnear, n] = find_nearpoint(Dist, cut, numb, chois)
% find n elements nearest to fixed point

b1 = -1;
b2 = -1;
ds = 0.1*cut;
dm = 0.01*cut;
if chois == 1
    for je = 1:1000
        cnear = find(Dist<=cut);
        n = length(cnear);
        if n > numb
            cut = cut*(1-ds);
            b1 = je;
        elseif n < numb
            cut = cut*(1+ds);
            b2 = je;
        end
        if b2-b1 == 1
            if ds >= dm+dm
                ds = ds-dm;
            else
                ds = ds*0.9;
            end
        end
        if n==numb || (ds<=0 && abs(b1-b2)==1)
            break;
        end
    end
    
else
    for je = 1:1000
        cnear = find(Dist>=cut);
        n = length(cnear);
        if n > numb
            cut = cut*(1+ds);
            b1 = je;
        elseif n < numb
            cut = cut*(1-ds);
            b2 = je;
        end
        if b2-b1 == 1
            if ds >= dm+dm
                ds = ds-dm;
            else
                ds = ds*0.9;
            end
        end
        if n==numb || (ds<=0 && abs(b1-b2)==1)
            break;
        end
    end
    
end
end

function [clusters, newlabels] = ind2cluster(labels)
% Transform vector of integers (cluster labels) to cell array.
C = unique(labels);
newlabels = labels;
k = length(C);
clusters = cell(1,k);
for i = 1:k
    ind = find(labels==C(i));
    clusters{i} = ind;
    newlabels(ind) = i;
end

end
