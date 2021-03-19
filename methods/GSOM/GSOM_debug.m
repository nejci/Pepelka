
function [retVal]=gSOM_debug(data,G,dG,alfa,maxIter,pOut,msize,shape,showSOM,showProcess,isMass,massThr,nClustPredef)

%   data            %podatki (morajo biti normalizirani na [0,1])
%   G=0.0008;       %vpliv gravitacije
%   dG=0.045;       %za koliko % se zmanjša G vsako iteracijo
%   alfa=0.01       %blizina, potrebna za zdruzitev
%   maxIter=100;    %stevilo iteracij zanke
%   pOut=0.1        %verjetnost, da namesto soseda izberemo nakljucni
%                   vzorec izmed preostalih, ko premikamo dani delec
%   msize=[]        %velikost SOM mreže, èe je [] se doloèi avtomatsko
%   shape='rect'    %oblika SOM mreze (rect ali hexa)
%   showSOM         %0 - ne prikazi SOM mreze na koncu 1. faze
%                   %1 - prikazi mrezo SOM
%   showProcess     %0 - ne prikazuj rezultatov
%                   %1 - prikazi vse vmesne korake
%                   %2 - prikazi koncni rezultat
%
%   isMass          %0 - ne upostevamo mase delcev
%                   %1 - upostevamo maso delcev - vecja masa, vecji privlak
%                   %2 - obratna gravitacija - vecja masa, manjsi privlak
%   massThr         %gledamo maso BMUjev (stevilo vzorcev, ki jih obsegajo)
%                   %in vse, ki imajo maso manjso od massThr, odstranimo
%
%   nClustPredef    %najmanjse zahtevano stevilo gruc
%
%     data=[  1,1;
%             2,1;
%             1,2;
%             2,2;
%             -2,1;
%             -2,2];
%     data=sintData(10,2,[-1,-1,1,1],0.6,0);
%
%Default:
%[retVal]=universeSuns(data,0.0008,0.045,0.01,100,0.1,[],'rect',1,2,0,2);

    %preusmerimo pot na spremenjeni som toolbox (hitrejsi PCA, popravek za rect/hexa size)
    addpath('..\..\misc\somtoolbox_pplk');
    

    eps=1e-3;
    %neighWeight=2.0;
    %verjetnost, da ni izbran eden od sosedov
    %pOut=0.1;
    
	%seznam iteracij, ki jih zelimo prikazati
	showIterList=[1,2,3,5,10,20,35];
	
    data=som_normalize(data, 'range');
    dataStart=data;
    
    [nPstart,dimStart]=size(dataStart);
    
    tic;
    
    [suns,sLabels]=createSOM(data,msize,shape,showSOM, isMass, massThr); 
    
      
    tSOM=toc;
    
    data=suns.coords;
    
    neighC=suns.neigh;
    mass=suns.mass;
    
    len=size(data,1);
    
    
    if(showProcess>0)
        f=figure();
        
        if dimStart > 2
            %izracunati moramo PCA
            if isempty(suns.dataPca)
                eig_vec=princomp(dataStart,'econ');
                dataStartOK=(eig_vec(:,1:2)' * dataStart')';
            else
                dataStartOK=suns.dataPca;
            end

             %pripravimo podatke za izris sonc vnaprej
            if isempty(suns.pca_vec)
                %uporabimo razultat PCA nad data
                suns_pca_vec=eig_vec; %princomp(suns.coords,'econ');
            else
                suns_pca_vec=suns.pca_vec;
            end
        else
            dataStartOK=dataStart;
        end
        %plot(dataStartOK(:,1),dataStartOK(:,2),'bo');
        axis('square')
        
        
       
        
    end
    %pause()
    %napravimo podatkovno strukturo za hranjenje strukture združevanja
    %imamo celièno strukturo -> dolgo len
    %ko združimo dva elementa, ju damo v isto celico in prazno izbrišemo
    %združena elementa se obnašata kot en vzorec (da se iznièi vpliv veèje gostote vzorcev)
    %pozicija združenega elementa je na poziciji prvega vzorca (?!?!?!?)
    
    %for c=1:len
        %C{c}={c,0}; %drugi element je 0, da se lahko na sklicujemo na prvi element kot na C{i}{1}
    %end
    %C je seznam vseh tock, ki se dopolnjuje pri zdruzevanju
    C=num2cell(1:len);
    
    %valid je tabela indeksov na tocke, ki so se veljavne
    valid=1:len;
    
    %-----------------------------------------------------
    %GLAVNA ZANKA
    nClust=len;
    
    %zastavica, ki pove, ali smo alg. koncali zaradi dosega koncnega stevila gruc
    retVal.nClustReached=0;
    
		
    tStart=tic;
    %frame=1;
    for iter=1:maxIter
		
		%zastavica z aprikazovanje samo prve slike v iteraciji
		showIter=1;
        
        %len=length(C); %dolzina se dinamicno spreminja (manjsa)
        i=1; %trenutno stanje stevca po strukturi C
        
        if nClust==1
            break;
        end
        
        %ali je stevilo gruc ze ustrezno zeljenemu stevilu?
        if nClust<=nClustPredef
            if showProcess
                fprintf('Exiting - nClustPredef reached!');
            end
            retVal.nClustReached=1;
            break;
        end
        
        %zastavica, ki signalizira konec premikanja delcev
        converged=1;
        
        if(showProcess > 0)
           disp(['Iteracija: ',num2str(iter),', G=',num2str(G),' ...']);
        end
        
        %napravimo nakljuèni vrstni red obiskovanja prototipov
               
        %while i <= len
		visitOrder=valid(randperm(length(valid)));
        for i=visitOrder
            %ce je i-ti element izbrisan, ga preskocimo
            if ~isempty(C{i})
                %izberemo nakljucni vzorec, ki ni i
                %prednost imajo sosedje, ce obstajajo

                %seznam sosedov i-tega sonca
                ni=neighC{i};

                if (~isempty(ni)) && (rand < 1-pOut)
                    %nakjlucno izberemo enega od sosedov
                    rpi=randperm(length(ni));
                    k=ni(rpi(1));

                else
                    %pozor! preveri, ko ne brisemo
                    %izberemo nakljucnega, lahko tudi soseda
                    choices=valid;
                    choices(valid==i)=[];
                    rpi=randperm(nClust-1);
                    k=choices(rpi(1));
                end


                %indeksi v tabeli data, ki so dostopni preko tabele C
                %dataI = C{i}{1};
                %dataK = C{k}{1};

                %vektor premika
                R=data(k,:)-data(i,:);
                E=R./norm(R);
                
                %nova pozicija i-tega in k-tega vzorca
                %TODO:premik lahko obtezimo glede na to, ali sta dva soseda ali ne
                
                %racunamo silo med dvema delcema, upostevamo maso ali pa ne
                if(isMass==1)
                    F=(R.*G*mass(i)*mass(k))./(norm(R)^3);
                %obratna gravitacija - vecja masa, manjsi privlak -
                %eksperimentalno
                elseif(isMass==2)
                    F=(R.*G)./(norm(R)^3 *mass(i)*mass(k));
                else
                    F=(R.*G)./(norm(R)^3);
                end
                %ce je velikost premika vecja od polovice razdalje, jo
                %omejimo s to mejo
                premik=min(norm(R)/2 , norm(F) );


                %varovalka za iteracije - ce ni bilo premikov, koncamo
                if premik > eps
                    converged=0;
                end

                data(i,:)=data(i,:)+ premik.*E;
                data(k,:)=data(k,:) - premik.*E;            

                %ali sta dovolj blizu skupaj za zdruzitev?
                if(norm(data(i,:)-data(k,:)) < alfa)
                    %k-ti element (z vsemi clani) pridružimo i-temu
                    C{i} = [C{i} C{k}];
                    %maso i-tega delca povecamo za maso k-tega
                    mass(i)=mass(i)+mass(k);
    
                    %morebitne 0 izbrišemo
                    %Carr=cell2mat(C{i});
                    %Carr(Carr==0)=[];
                    %C{i}=num2cell(Carr);

                    %k-ti element postavimo na [], kar pomeni neveljaven
                    C{k}=[];
                    valid(valid==k)=[];
                    %vodimo evidenco o preostalih soncih
                    nClust=nClust-1;

                    %-----------------------------------------------
                    %popraviti moramo podatke o sosedih
                    %k odstranimo iz seznama sosedov i-ja
                    new=neighC{i};
                    new(new==k)=[];                
                    neighC{i}=new;

                    %vsem sosedom k-ja prevezemo reference na i
                    %seznam sosedov k-ja
                    nk=neighC{k};

                    for nInd=1:length(nk)
                       new=neighC{nk(nInd)};
                       if ~ismember(i,new)
                         new(new==k)=i;
                       else
                         new(new==k)=[]; 
                       end
                       neighC{nk(nInd)}=new;
                       
                       %sosedi k-ja postanejo tudi sosedi i-ja
                       if ~ismember(nk(nInd),neighC{i}) && nk(nInd) ~=i
                           neighC{i}=[neighC{i}, nk(nInd)];
                       end
                    end

                    %k damo na [], kar pomeni,da nima sosedov
                    neighC{k}=[];

                end

                %fprintf('Par i=%d, k=%d\n',i,k);

                if(showProcess==1 && showIter==1 && ismember(iter,showIterList))
					
					%pomeni, da bo prikazana samo ena slika na zacetku
					%iteracije in ne premik vsake tocke posebej.
					showIter=0;
					
					figure();
					
                    %pause(0.005);
                    %pause();
                    %plot(dataStartOK(:,1),dataStartOK(:,2),'bo');
                    axis('square')
                    hold on;
                    
                    %na podlagi prvih dveh lastnih vektorjev izracunamo
                    %nove pozicije sonc
                    if dimStart > 2
                        dataOK=(suns_pca_vec(:,1:2)' * data')';
                        
                    else
                        dataOK=data;
                    end
                    
                    for j=1:nClust
                        %if ~isempty(C{j})
                            plot(dataOK(C{valid(j)}(1),1),dataOK(C{valid(j)}(1),2),'k.','markersize',10);
                            %plot(dataOK(C{valid(j)}(1),1),dataOK(C{valid(j)}(1),2),'r.');
                            
							for sosedI=1:length(neighC{C{valid(j)}(1)})
								izvorIData= C{valid(j)}(1);
								sosedIData=neighC{C{valid(j)}(1)}(sosedI);
								
								if izvorIData < sosedIData 
									plot(	[	dataOK(izvorIData,1); dataOK(sosedIData,1)], ...
											[	dataOK(izvorIData,2);dataOK(sosedIData,2)],'k-');
							
								end
							end
							
							title(['Iteracija: ', num2str(iter),', G=',num2str(G)]);
							
							%text(dataOK(C{valid(j)}(1),1),dataOK(C{valid(j)}(1),2),num2str(mass(C{valid(j)}(1))),'verticalAlignment','bottom');
                            %M(frame)=getframe;
                            %frame=frame+1;
                        %end
                    end
                    hold off;
                    if dimStart <=2
                        axis([0,1,0,1]);
                    end
                    
                end
                
            end
            %i=i+1;
        end
        %ce je bil vsak premik pod pragom eps
        if converged
            break;
        end
        G=(1-dG)*G;
        
    end
    %GLAVNA ZANKA - END
    %save('movie.mat','M');
    %-----------------------------------------------------
    
    %pause()
    
    %----------------------------------------------------
    %rezultati clusteringa
    %nClust=length(C);
    target=zeros(size(dataStart,1),1);
    
      
    %restavracija oznak za originalne podatke
    for c=1:nClust
        
        tmp=C{valid(c)};
        %tmp(tmp==0)=[];

        sL=suns.labels(tmp);

        for s=1:length(sL)
            target(sL(s)==sLabels)=c;
        end
        
    end
    
    tG=toc;
    
    if(showProcess==2)
        plot(dataStartOK(:,1),dataStartOK(:,2),'ro');
        hold on;
        
        if dimStart > 2
            dataOK=(suns_pca_vec(:,1:2)' * data')';
        else
            dataOK=data;
        end
        
        for j=1:nClust
            plot(dataOK(C{valid(j)}(1),1),dataOK(C{valid(j)}(1),2),'.');
        end
        hold off;
        if dimStart <=2
            axis([0,1,0,1]);
        end
                    
    end
    
    data=dataStart;
    
    if(showProcess>0)
       draw_clustering2(data,nClust,target,'gSOM',true,1);        
    end
    
    retVal.target=target;
    retVal.nClust=nClust;
    retVal.iter=iter;
    retVal.time=[tSOM, tG, tSOM+tG];
    retVal.SOMprop=suns.SOMprop;
    retVal.SOMtrain=suns.SOMtrain;
    
    %izbrisemo pot na spremenjeni som toolbox
    rmpath('..\..\misc\somtoolbox_pplk');
    
end