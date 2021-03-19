% test of FEKM, FSKM

k = 4;

% load raw values of cluster validity indices
load CVI.mat

list_CVI_sel = {'DN','DB','SDBW','CH','SIL','I','XB'};
GDIpos = find(strcmp('GDI',list_CVI_sel));
if GDIpos
    GDInames = {'GDI_11','GDI_12','GDI_13',...
        'GDI_21','GDI_22','GDI_23',...
        'GDI_31','GDI_32','GDI_33',...
        'GDI_41','GDI_42','GDI_43',...
        'GDI_51','GDI_52','GDI_53',...
        'GDI_61','GDI_62','GDI_63'};
    CVI_header_sel = [list_CVI_sel(1:GDIpos-1),GDInames,list_CVI_sel(GDIpos+1:end)];
else
    CVI_header_sel = list_CVI_sel;
end

[Lia, Locb] = ismember(CVI_header, CVI_header_sel);
ind_CVI(Locb(Locb>0)) = find(Lia);
num_CVI = length(CVI_header_sel);


% take one matrix as sample
PRM_raw = CVI{1,1}(:,ind_CVI);



% unify
listOfMinLike = {'APN', 'AD', 'ADM', 'CI', 'CON', 'DB', 'DB*', 'FOM', ...
    'GAMMA','SD','SDBW','SEP','SEPMAX','TAU','VAR','XB','XB*'};
mask = ismember(CVI_header_sel,listOfMinLike);

PRM = pplk_unifyPRM(PRM_raw,mask,'minmax');

% compute FEKM
[PRM_FEKM,labels_FEKM,energy_FEKM] = FEKM(PRM,k);
fprintf(1,'FEKM energy: %f\n',energy_FEKM);

% compute FSKM
[PRM_FSKM,featInd_FSKM,labels_FSKM,energy_FSKM] = FSKM(PRM,k);
fprintf(1,'FSKM energy: %f\n',energy_FSKM);
fprintf(1,'FSKM selects: ');fprintf(1,'%s; ', CVI_header_sel{featInd_FSKM}); fprintf(1,'\n');

% plot FEKM
figure();
plot(PRM,'k-.'); hold on; plot(PRM_FEKM,'r'); title('FEKM');
figure();
pplk_scatterPlot(PRM',labels_FEKM);

% plot FEKM
figure();
plot(PRM,'k-.'); hold on; plot(PRM_FSKM,'r');title('FSKM');
figure();
pplk_scatterPlot(PRM',labels_FSKM);


% compare with LS, FSFS and SPEC
[PRM_LS,featInd_LS] = pplk_featureReduce(PRM, 'LS', k);
fprintf(1,'LS selects: ');fprintf(1,'%s; ', CVI_header_sel{featInd_LS}); fprintf(1,'\n');
figure(); plot(PRM,'k-.'); hold on; plot(PRM_LS,'r'); title('LS');


[PRM_FSFS,featInd_FSFS] = pplk_featureReduce(PRM, 'FSFS', k);
fprintf(1,'FSFS selects: ');fprintf(1,'%s; ', CVI_header_sel{featInd_FSFS}); fprintf(1,'\n');
figure(); plot(PRM,'k-.'); hold on; plot(PRM_FSFS,'r'); title('FSFS');

[PRM_SPEC,featInd_SPEC] = pplk_featureReduce(PRM, 'SPEC', k);
fprintf(1,'SPEC selects: ');fprintf(1,'%s; ', CVI_header_sel{featInd_SPEC}); fprintf(1,'\n');
figure(); plot(PRM,'k-.'); hold on; plot(PRM_SPEC,'r'); title('SPEC');



