
function [suns,labels]=createSOM(data,msize,shape,show,isMass,massThr)
%[suns,labels]=createSuns(data)
%
%Phase I. - aproximating data with SOM
%
%INPUTS
%   data : d-dimensional data on which SOM is trained
%   
%
%OUTPUTS:
%   suns    :   stucture with fields
%                   .labels : labels of SOM units [to link with data patterns]
%                   .coords : coordinates of SOM units
%                   .neigh  :   pairs of neighbours (a,b), where a < b
%                   .mass : mass (number of particles) of each BMU
%
%   labels  :   labels of BMUs - each data pattern has its
%               coresponding BMU.
%   


[nP,dim]=size(data);

D=som_data_struct(data);

%normalizacija na [0,1] - daj ven, ker to ze predpostavljamo
%D = som_normalize(D, 'range');

%som_make privzame batch nacin ucenja
%sMap = som_make(D, 'msize', [3,3], 'lattice', 'hexa','shape', 'sheet'); 

if strcmp(shape,'hexa')
    if isempty(msize)
        sMap = som_make(D,'lattice', 'hexa','shape', 'sheet', 'tracking',0);         
    else
        sMap = som_make(D,'msize',msize,'lattice', 'hexa','shape', 'sheet', 'tracking',0); 
    end
else
    if isempty(msize)
        sMap = som_make(D, 'lattice', 'rect','shape', 'sheet', 'tracking',0); 
    else
        sMap = som_make(D, 'msize',msize,'lattice', 'rect','shape', 'sheet', 'tracking',0); 
    end
end
%sMap = som_make(D); 

%sMap=som_autolabel(sMap, D, 'vote');

suns.SOMprop=sMap.topol;
suns.SOMtrain=sMap.trainhist;

%Compute BMUs (1st, 2nd and 3rd) for data patterns
%v primeru, da zaradi mase dolocen BMU izbrisemo, ga lahko hitro zamenjamo
%z drugim ali tretjim
which=1:3;
labels=som_bmus(sMap, D.data, which);
labelsOth=labels(:,2:length(which));
labels=labels(:,1);

%determine unique BMUs -> these are our suns in the galaxy
labelsUniq=unique(labels);

%how many data points are represented by certain BMU? BMU labels are
%sorted.
mass=hist(labels,labelsUniq);


%[labels,labelsUniq,mass]=massThresholdCut(data,labels,labelsOth,labelsUniq,mass,massThr,sMap);
%[labels,labelsUniq,mass]=massThresholdIncr(data,labels,labelsOth,labelsUniq,mass,massThr,sMap);

%shranimo podatke v strukturo suns
suns.labels=labelsUniq;
suns.coords = sMap.codebook(labelsUniq,:);

if(isMass)
    suns.mass=mass;
else
    suns.mass=ones(1,length(mass));
end


%compute neighbours
N=som_connection(sMap);


%create matrix Npairs x 2 for storing neighbour info
pairs=nchoose2(labelsUniq);

neigh=[];

%check whether certain pair is neighbouring
for ind=1:size(pairs,1)
    if(N(pairs(ind,1),pairs(ind,2))==1)
        neigh=[neigh; pairs(ind,:)];
    end
    
end

%create cell struct - each sun has a list of its neighbours (indices in labelsUniq)
neighC=cell(1,length(labelsUniq));
for s=1:length(labelsUniq)
    if isempty(neigh)
        x=[];
        y=[];
    else
        x=find(neigh(:,1)==labelsUniq(s));
        y=find(neigh(:,2)==labelsUniq(s));
    end
    
    if ~isempty(x)
        nX=neigh(x,2);
        for nXi=1:length(nX)
            neighC{s}=[neighC{s}, find(labelsUniq==nX(nXi))];
        end
    end
    if ~isempty(y)
        nY=neigh(y,1);
        for nYi=1:length(nY)
            neighC{s}=[neighC{s}, find(labelsUniq==nY(nYi))];
        end
    end
end

suns.neigh=neighC;

%colormap(1-gray);
%som_show(sMap, 'norm', 'd');

%lastni vektorji sonc - to rabimo pozneje za izris
suns.pca_vec=[];
%ze preslikani podatki na 2 dimenziji
suns.dataPca=[];

if(show)
        
    if dim > 2
        %PCA preslikava nevronov in sonc
        neur=[sMap.codebook ; suns.coords];
        %eig_vec_CB=princomp(neur,'econ');
        %new_CB=(eig_vec_CB(:,1:2)' * neur')';
        %kot glavni komponenti uporabimo kar tisti, izracunani iz podatkov
        %data
        new_CB=(D.pca_vec(:,1:2)' * neur')';
        
        
        numCB=size(sMap.codebook,1);

        CB=new_CB(1:numCB,:);
        suns_new=new_CB(numCB+1:end,:);
        
        %PCA preslikava podatkov
        new_data=(D.pca_vec(:,1:2)' * D.data')';
        
        suns.pca_vec=D.pca_vec; %eig_vec_CB;
        suns.dataPca=new_data;
        
    else
        CB=sMap.codebook;
        suns_new=suns.coords;
        new_data=D.data;
        
    end
    
    figure();
    som_grid(sMap, 'Coord', [CB(:,1),CB(:,2)]);
    axis('equal')
    hold on;
    
    plot(new_data(:,1),new_data(:,2),'bo');
    
    for i=1:length(suns.mass)
        plot(suns_new(i,1),suns_new(i,2),'r.', 'markersize', suns.mass(i)/max(suns.mass)*20)
        %text(suns_new(i,1),suns_new(i,2)+0.01,num2str(suns.mass(i)),'verticalAlignment','bottom');
    end
    hold off;
    
end
end %function

function [labels,labelsUniq,mass]=massThresholdIncr(data,labels,labelsOth,labelsUniq,mass,massThr,sMap)
% odstranitev nevronov glede na prag mase - inkrementalno. Najprej za
% prag=1, 2, ... , massThr
    
    %prag za maso - nevron je sprejet, ce ima maso vecjo ali enako tej
    %kako doloèiti ta prag?
    
    for massThr_i=1:massThr

    
        %seznam vseh nevronov BMU, ki jih je potrebno odstraniti, ker imajo
        %premajno maso
        blackList=labelsUniq(mass < massThr_i); %%oznake nevronov v seznamu labelsUniq
        blackListI=find(mass < massThr_i); %indeksi nevronov v seznamu labelsUniq

        %seznam indeksov sirot - to so tisti vzorci, ki so imeli za BMU nekoga iz
        %seznama blackList
        orphans=[];
        for ind=1:length(blackList)
            orphans=[orphans; find(labels==blackList(ind))];
        end

        %odstranimo BMUje s premajhno maso
        labelsUniq(blackListI)=[];
        mass(blackListI)=[];

        coords = sMap.codebook(labelsUniq,:);


        %popravimo labels za vzorce - sirotam poiscemo starse
        %ce ima sirota v labelsOth nevron, ki ne odpade, se ta uporabi, sicer se
        %poisce najblizjega.
        %zanka po vseh sirotah
        for orph=1:length(orphans)
            %zanka po vseh sekundarnih BMUjih 
            %ali smo nasli nove starse?
            parentFound=false;

            for b=1:size(labelsOth,2)
                %ce najdemo tak sekundarni BMU, ki je eden od preostalih z dovolj veliko maso,
                %ga sprejmemo kot novi BMU za siroto
                answer=(labelsOth(orphans(orph),b) == labelsUniq);
                if sum(answer)==1
                    labels(orphans(orph))=labelsOth(orphans(orph),b);
                    %masa se poveca za 1!
                    mass(answer)=mass(answer)+1;

                    parentFound=true;
                    break;
                end       
            end

            %ce noben sekundarni BMU ni ustrezen, ga poiscemo sami
            %izracunamo evklidsko razdaljo sirote do vseh nevronov, ki so v labelsUniq
            if ~parentFound        
                disp('iscem starse-razdalja')
                a=coords;
                b=data(orphans(orph),:);
                d=(sum( (a-repmat(b,size(a,1),1)).^2,2)).^0.5;

                [v,ind]=min(d);

                labels(orphans(orph))=labelsUniq(ind);
                %masa se poveca za 1!
                mass(ind)=mass(ind)+1;
            end
        end
        
        
        %izracun cenilke - TODO
    end
    
end

function [labels,labelsUniq,mass]=massThresholdCut(data,labels,labelsOth,labelsUniq,mass,massThr,sMap)
% odstranitev nevronov glede na prag mase - hkratno. Takoj za dolocen prag.
    
    %te spremenljivke se spremenijo, zapomnimo si originale
    labelsUniqOrig=labelsUniq;
    labelsOrig=labels;
    massOrig=mass;
    
    internal_history=[];
    wtertra_h=0; %zacetna vrednost, da ni problema na koncu pri sklicevanju na indeks end-1
    
    for massThr_i=1:massThr
    
        %shranimo podatke od prejsnjega kroga
        labelsUniq_old=labelsUniq;
        labels_old=labels;
        mass_old=mass;
        
        %resetiramo podatke na originalne vrednosti
        labelsUniq=labelsUniqOrig;
        labels=labelsOrig;
        mass=massOrig;
        
        %seznam vseh nevronov BMU, ki jih je potrebno odstraniti, ker imajo
        %premajno maso
        blackList=labelsUniq(mass < massThr_i); %%oznake nevronov v seznamu labelsUniq
        blackListI=find(mass < massThr_i); %indeksi nevronov v seznamu labelsUniq

        %seznam indeksov sirot - to so tisti vzorci, ki so imeli za BMU nekoga iz
        %seznama blackList
        orphans=[];
        for ind=1:length(blackList)
            orphans=[orphans; find(labels==blackList(ind))];
        end

        %odstranimo BMUje s premajhno maso
        labelsUniq(blackListI)=[];
        mass(blackListI)=[];

        coords = sMap.codebook(labelsUniq,:);


        %popravimo labels za vzorce - sirotam poiscemo starse
        %ce ima sirota v labelsOth nevron, ki ne odpade, se ta uporabi, sicer se
        %poisce najblizjega.
        %zanka po vseh sirotah
        for orph=1:length(orphans)
            %zanka po vseh sekundarnih BMUjih 
            %ali smo nasli nove starse?
            parentFound=false;

            for b=1:size(labelsOth,2)
                %ce najdemo tak sekundarni BMU, ki je eden od preostalih z dovolj veliko maso,
                %ga sprejmemo kot novi BMU za siroto
                answer=(labelsOth(orphans(orph),b) == labelsUniq);
                if sum(answer)==1
                    labels(orphans(orph))=labelsOth(orphans(orph),b);
                    %masa se poveca za 1!
                    mass(answer)=mass(answer)+1;

                    parentFound=true;
                    break;
                end       
            end

            %ce noben sekundarni BMU ni ustrezen, ga poiscemo sami
            %izracunamo evklidsko razdaljo sirote do vseh nevronov, ki so v labelsUniq
            if ~parentFound        
                disp('iscem starse-razdalja')
                a=coords;
                b=data(orphans(orph),:);
                d=(sum( (a-repmat(b,size(a,1),1)).^2,2)).^0.5;

                [v,ind]=min(d);

                labels(orphans(orph))=labelsUniq(ind);
                %masa se poveca za 1!
                mass(ind)=mass(ind)+1;
            end
        end
        
        
           
        %izracunamo cenilko
        retVal=valid_internal(data,labels,1);
        internal_history=[internal_history retVal];
        wtertra_h=[wtertra_h retVal.Wtertra];
        wtertra_h(isnan(wtertra_h))=0; %morebitne NaN vrednosti damo na 0
        %ce je indeks Wtertra za trenutni massThr_i manjsi od prejsnjega,
        %so vrednosti prejsnje iteracije tiste optimalne. ALI
        %v primeru, da ni sprememb od prejsnje iteracije (pomeni, da ni
        %nobenih novih za odstraniti), koncamo ponavljanje
        
        if (wtertra_h(end-1)>wtertra_h(end)) %|| isempty(blackList))
           disp('Nasli smo ga!');
           
           %vrnemo vrednosti pri optimalnem pragu
           labelsUniq=labelsUniq_old;
           labels=labels_old;
           mass=mass_old;
           
           break;
        end
        
    end
   
end

