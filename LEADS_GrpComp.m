%%% Script to automatically generate group comparisons and
%%% clinico-functional correlations using as much as data possible in LEADS
%%% Leonardo.Iaccarino@ucsf.edu October 2020 - LEADS.PETCORE@ucsf.edu

% set paths that could be useful
addpath('/mnt/coredata/Projects/LEADS/script_f7p1/');
addpath('/mnt/coredata/Projects/LEADS/script_f7p1/service/');
addpath('/mnt/coredata/Projects/LEADS/script_f7p1/templates/');
addpath('/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/');
addpath(genpath('/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/'));
addpath(genpath('/mnt/neuroimaging/SPM/spm12/matlabbatch/'));
path_master='/mnt/coredata/Projects/LEADS/data_f7p1/';
path_compcor='/mnt/coredata/Projects/LEADS/data_f7p1/stats/';
path_compcor_grpcomp='/mnt/coredata/Projects/LEADS/data_f7p1/stats/grpcomp/';
path_dbs='/shared/petcore/Projects/LEADS/data_f7p1/dbs/';
path_clindf='/shared/petcore/Projects/LEADS/data_f7p1/LONI_uploads/ready_datasets/';

emask='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring/mask/EM.nii';

% need to read latest database with both clinical and imaging data 

dbclinim=dir(strcat(path_clindf,'FBBFTP_Baseline*')); % dir for the FBBFTP databases
dbclinim=struct2cell(dbclinim)';
dbclinim=dbclinim(:,[1 3]);
dbclinim=array2table(dbclinim);
dbclinim.Properties.VariableNames([1]) = cellstr(strcat('Filename'));
dbclinim.Properties.VariableNames([2]) = cellstr(strcat('DateCreated'));
dbclinim.DateCreated=datetime(dbclinim.DateCreated);
filt=max(dbclinim.DateCreated); % Store most recent date
dbclinim=dbclinim(dbclinim.DateCreated==filt,:); % Subset to select the most recent file
dbclinimf=readtable(strcat(path_clindf,char(table2cell(dbclinim(1,1))))); %% loading the clin database

% need to read latest path database

dbimpaths=dir(strcat(path_dbs,'LEADS_DatabasePaths*')); % dir for the FBBFTP databases
dbimpaths=struct2cell(dbimpaths)';
dbimpaths=dbimpaths(:,[1 3]);
dbimpaths=array2table(dbimpaths);
dbimpaths.Properties.VariableNames([1]) = cellstr(strcat('Filename'));
dbimpaths.Properties.VariableNames([2]) = cellstr(strcat('DateCreated'));
dbimpaths.DateCreated=datetime(dbimpaths.DateCreated);
filt=max(dbimpaths.DateCreated); % Store most recent date
dbimpaths=dbimpaths(dbimpaths.DateCreated==filt,:); % Subset to select the most recent file
dbimpathsf=readtable(strcat(path_dbs,char(table2cell(dbimpaths(1,1))))); %% loading the clin database

% now merge the two databases 

dbm=innerjoin(dbclinimf,dbimpathsf);

% Remove the two Amy- AD-Tau+, LDS1770072 and LDS0370175
% Update December 7th, now three cases considering LDS0360311
% Update April 9th, now four cases including LDS0360322 very positive CN

dbm(strcmp(dbm.ID,'LDS1770072'),:)=[];
dbm(strcmp(dbm.ID,'LDS0370175'),:)=[];
dbm(strcmp(dbm.ID,'LDS0360311'),:)=[];
dbm(strcmp(dbm.ID,'LDS0360322'),:)=[];

% Time to organize group comparisons

dir_grpcomp=strcat(path_compcor_grpcomp,datestr(now,'mm-dd-yyyy'));
mkdir (dir_grpcomp);

% We need to subset the merged table to organize the pairwise comparisons

dbm_eoadeononad=dbm(~strcmp(dbm.CohortAssgn_x, 'CN'), :);
dbm_eoadcn=dbm(~strcmp(dbm.CohortAssgn_x, 'EOnonAD'), :);
dbm_eononadcn=dbm(~strcmp(dbm.CohortAssgn_x, 'EOAD'), :);
dbm_eoad=dbm(strcmp(dbm.CohortAssgn_x, 'EOAD'), :); % we can use this later for EOAD APOE e4+ vs. APOE e4-

% Select variables we are interested in

dbm_eoadeononad=dbm_eoadeononad(:,[1 9 2 14 35 24 26 28 30 42 46]);
dbm_eoadcn=dbm_eoadcn(:,[1 9 2 14 35 24 26 28 30 42 46]);
dbm_eononadcn=dbm_eononadcn(:,[1 9 2 14 35 24 26 28 30 42 46]);
dbm_eoad=dbm_eoad(:,[1 9 2 14 35 24 26 28 30 42 46]);

% cleanup databases from NaNs

dbm_eoadeononad=dbm_eoadeononad(~any(ismissing(dbm_eoadeononad),2),:);
dbm_eoadcn=dbm_eoadcn(~any(ismissing(dbm_eoadcn),2),:);
dbm_eononadcn=dbm_eononadcn(~any(ismissing(dbm_eononadcn),2),:);
dbm_eoad=dbm_eoad(~any(ismissing(dbm_eoad),2),:);

% recode sex variable in a dummy

dbm_eoadeononad.sex_rcd=strcmp(dbm_eoadeononad.sex,'Female');
dbm_eoadcn.sex_rcd=strcmp(dbm_eoadcn.sex,'Female');
dbm_eononadcn.sex_rcd=strcmp(dbm_eononadcn.sex,'Female');
dbm_eoad.sex_rcd=strcmp(dbm_eoad.sex,'Female');

% sort by group

dbm_eoadeononad=sortrows(dbm_eoadeononad,{'CohortAssgn_x'},{'ascend'});
dbm_eoadcn=sortrows(dbm_eoadcn,{'CohortAssgn_x'},{'descend'});
dbm_eononadcn=sortrows(dbm_eononadcn,{'CohortAssgn_x'},{'descend'});
dbm_eoad=sortrows(dbm_eoad,{'APOEpos'},{'descend'});

% ready to build the SPM group comparisons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% 1st: EOAD vs. EOnonAD %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_grpcomp_eoadeononad=strcat(dir_grpcomp,'/EOADvsEOnonAD_agesexmmseapoe'); mkdir (dir_grpcomp_eoadeononad);
dir_grpcomp_eoadeononad_fbb=strcat(dir_grpcomp,'/EOADvsEOnonAD_agesexmmseapoe/fbb'); mkdir (dir_grpcomp_eoadeononad_fbb);
dir_grpcomp_eoadeononad_ftp=strcat(dir_grpcomp,'/EOADvsEOnonAD_agesexmmseapoe/ftp'); mkdir (dir_grpcomp_eoadeononad_ftp);

n_eoad=nnz(strcmp(dbm_eoadeononad.CohortAssgn_x,'EOAD'));
n_eononad=size(dbm_eoadeononad,1)-n_eoad;

% Print a small log
            
Vars={'N EOAD'; 'N EOnonAD'; 'Covariates'};
Vals={n_eoad; n_eononad; 'Age,Sex,MMSE,APOE'};
tempT=table(Vars,Vals);
writetable(tempT, strcat(dir_grpcomp_eoadeononad,'/Log.csv'));

% FBB

eoad_wfbb=table2cell(dbm_eoadeononad(1:n_eoad,{'Timepoint1_Path_FBB_wsuvr'}));
eononad_wfbb=table2cell(dbm_eoadeononad(n_eoad+1:size(dbm_eoadeononad,1),{'Timepoint1_Path_FBB_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eoadeononad_fbb);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = eoad_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = eononad_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eoadeononad.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eoadeononad.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c =  dbm_eoadeononad.MMSE;
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'mmse';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).c = str2double(dbm_eoadeononad.APOEpos);
matlabbatch{1}.spm.stats.factorial_design.cov(4).cname = 'apoe';
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'EOAD>EOnonAD';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'EOAD<EOnonAD';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'eoadMeononad_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'eoadMeononad_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'eoadLeononad_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'eoadLeononad_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_fbb,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_fbb,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_fbb,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_fbb,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eoadeononad_fbb,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eoadeononad_fbb,'/D_EOADgtEOnonAD_pFWE*')); es1p=strcat(dir_grpcomp_eoadeononad_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eoadeononad_fbb,'/D_EOADltEOnonAD_pFWE*')); es2p=strcat(dir_grpcomp_eoadeononad_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eoadeononad_fbb,'/D_EOADgtEOnonAD_p0.1*')); es1p=strcat(dir_grpcomp_eoadeononad_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eoadeononad_fbb,'/D_EOADltEOnonAD_p0.1*')); es2p=strcat(dir_grpcomp_eoadeononad_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end


clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eoadeononad)

% FTP

eoad_wftp=table2cell(dbm_eoadeononad(1:n_eoad,{'Timepoint1_Path_FTP_wsuvr'}));
eononad_wftp=table2cell(dbm_eoadeononad(n_eoad+1:size(dbm_eoadeononad,1),{'Timepoint1_Path_FTP_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eoadeononad_ftp);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = eoad_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = eononad_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eoadeononad.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eoadeononad.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c =  dbm_eoadeononad.MMSE;
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'mmse';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).c = str2double(dbm_eoadeononad.APOEpos);
matlabbatch{1}.spm.stats.factorial_design.cov(4).cname = 'apoe';
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'EOAD>EOnonAD';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'EOAD<EOnonAD';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'eoadMeononad_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'eoadMeononad_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'eoadLeononad_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'eoadLeononad_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_ftp,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_ftp,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_ftp,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_ftp,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eoadeononad_ftp,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eoadeononad_ftp,'/D_EOADgtEOnonAD_pFWE*')); es1p=strcat(dir_grpcomp_eoadeononad_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eoadeononad_ftp,'/D_EOADltEOnonAD_pFWE*')); es2p=strcat(dir_grpcomp_eoadeononad_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eoadeononad_ftp,'/D_EOADgtEOnonAD_p0.1*')); es1p=strcat(dir_grpcomp_eoadeononad_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eoadeononad_ftp,'/D_EOADltEOnonAD_p0.1*')); es2p=strcat(dir_grpcomp_eoadeononad_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% 2nd: EOAD vs. CN %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_grpcomp_eoadcn=strcat(dir_grpcomp,'/EOADvsCN_agesexmmseapoe'); mkdir (dir_grpcomp_eoadcn);
dir_grpcomp_eoadcn_fbb=strcat(dir_grpcomp,'/EOADvsCN_agesexmmseapoe/fbb'); mkdir (dir_grpcomp_eoadcn_fbb);
dir_grpcomp_eoadcn_ftp=strcat(dir_grpcomp,'/EOADvsCN_agesexmmseapoe/ftp'); mkdir (dir_grpcomp_eoadcn_ftp);

n_eoad=nnz(strcmp(dbm_eoadcn.CohortAssgn_x,'EOAD'));
n_cn=size(dbm_eoadcn,1)-n_eoad;

% Print a small log
            
Vars={'N EOAD'; 'N CN'; 'Covariates'};
Vals={n_eoad; n_cn; 'Age,Sex,MMSE,APOE'};
tempT=table(Vars,Vals);
writetable(tempT, strcat(dir_grpcomp_eoadcn,'/Log.csv'));

% FBB

eoad_wfbb=table2cell(dbm_eoadcn(1:n_eoad,{'Timepoint1_Path_FBB_wsuvr'}));
cn_wfbb=table2cell(dbm_eoadcn(n_eoad+1:size(dbm_eoadcn,1),{'Timepoint1_Path_FBB_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eoadcn_fbb);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = eoad_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = cn_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eoadcn.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eoadcn.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c =  dbm_eoadcn.MMSE;
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'mmse';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).c = str2double(dbm_eoadcn.APOEpos);
matlabbatch{1}.spm.stats.factorial_design.cov(4).cname = 'apoe';
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'EOAD>CN';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'EOAD<CN';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'eoadMcn_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'eoadMcn_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'eoadLcn_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'eoadLcn_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_fbb,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_fbb,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_fbb,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_fbb,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eoadcn_fbb,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eoadcn_fbb,'/D_EOADgtCN_pFWE*')); es1p=strcat(dir_grpcomp_eoadcn_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eoadcn_fbb,'/D_EOADltCN_pFWE*')); es2p=strcat(dir_grpcomp_eoadcn_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eoadcn_fbb,'/D_EOADgtCN_p0.1*')); es1p=strcat(dir_grpcomp_eoadcn_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eoadcn_fbb,'/D_EOADltCN_p0.1*')); es2p=strcat(dir_grpcomp_eoadcn_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eoadcn)

% FTP

eoad_wftp=table2cell(dbm_eoadcn(1:n_eoad,{'Timepoint1_Path_FTP_wsuvr'}));
cn_wftp=table2cell(dbm_eoadcn(n_eoad+1:size(dbm_eoadcn,1),{'Timepoint1_Path_FTP_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eoadcn_ftp);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = eoad_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = cn_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eoadcn.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eoadcn.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c =  dbm_eoadcn.MMSE;
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'mmse';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).c = str2double(dbm_eoadcn.APOEpos);
matlabbatch{1}.spm.stats.factorial_design.cov(4).cname = 'apoe';
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'EOAD>CN';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'EOAD<CN';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'eoadMcn_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'eoadMcn_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'eoadLcn_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'eoadLcn_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_ftp,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_ftp,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_ftp,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_ftp,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eoadcn_ftp,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eoadcn_ftp,'/D_EOADgtCN_pFWE*')); es1p=strcat(dir_grpcomp_eoadcn_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eoadcn_ftp,'/D_EOADltCN_pFWE*')); es2p=strcat(dir_grpcomp_eoadcn_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eoadcn_ftp,'/D_EOADgtCN_p0.1*')); es1p=strcat(dir_grpcomp_eoadcn_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eoadcn_ftp,'/D_EOADltCN_p0.1*')); es2p=strcat(dir_grpcomp_eoadcn_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eoadcn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% 3rd: EOnonAD vs. CN %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_grpcomp_eononadcn=strcat(dir_grpcomp,'/EOnonADvsCN_agesexmmseapoe'); mkdir (dir_grpcomp_eononadcn);
dir_grpcomp_eononadcn_fbb=strcat(dir_grpcomp,'/EOnonADvsCN_agesexmmseapoe/fbb'); mkdir (dir_grpcomp_eononadcn_fbb);
dir_grpcomp_eononadcn_ftp=strcat(dir_grpcomp,'/EOnonADvsCN_agesexmmseapoe/ftp'); mkdir (dir_grpcomp_eononadcn_ftp);

n_eononad=nnz(strcmp(dbm_eononadcn.CohortAssgn_x,'EOnonAD'));
n_cn=size(dbm_eononadcn,1)-n_eononad;

% Print a small log
            
Vars={'N EOnonAD'; 'N CN'; 'Covariates'};
Vals={n_eononad; n_cn; 'Age,Sex,MMSE,APOE'};
tempT=table(Vars,Vals);
writetable(tempT, strcat(dir_grpcomp_eononadcn,'/Log.csv'));

% FBB

eononad_wfbb=table2cell(dbm_eononadcn(1:n_eononad,{'Timepoint1_Path_FBB_wsuvr'}));
cn_wfbb=table2cell(dbm_eononadcn(n_eononad+1:size(dbm_eononadcn,1),{'Timepoint1_Path_FBB_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eononadcn_fbb);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = eononad_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = cn_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eononadcn.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eononadcn.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c =  dbm_eononadcn.MMSE;
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'mmse';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).c = str2double(dbm_eononadcn.APOEpos);
matlabbatch{1}.spm.stats.factorial_design.cov(4).cname = 'apoe';
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'EOnonAD>CN';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'EOnonAD<CN';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'eononadMcn_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'eononadMcn_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'eononadLcn_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'eononadLcn_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_fbb,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_fbb,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_fbb,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_fbb,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eononadcn_fbb,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eononadcn_fbb,'/D_EOnonADgtCN_pFWE*')); es1p=strcat(dir_grpcomp_eononadcn_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eononadcn_fbb,'/D_EOnonADltCN_pFWE*')); es2p=strcat(dir_grpcomp_eononadcn_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eononadcn_fbb,'/D_EOnonADgtCN_p0.1*')); es1p=strcat(dir_grpcomp_eononadcn_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eononadcn_fbb,'/D_EOnonADltCN_p0.1*')); es2p=strcat(dir_grpcomp_eononadcn_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eononadcn)

% FTP

eononad_wftp=table2cell(dbm_eononadcn(1:n_eononad,{'Timepoint1_Path_FTP_wsuvr'}));
cn_wftp=table2cell(dbm_eononadcn(n_eononad+1:size(dbm_eononadcn,1),{'Timepoint1_Path_FTP_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eononadcn_ftp);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = eononad_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = cn_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eononadcn.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eononadcn.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c =  dbm_eononadcn.MMSE;
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'mmse';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).c = str2double(dbm_eononadcn.APOEpos);
matlabbatch{1}.spm.stats.factorial_design.cov(4).cname = 'apoe';
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'EOnonAD>CN';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'EOnonAD<CN';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'eononadMcn_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'eononadMcn_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'eononadLcn_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'eononadLcn_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_ftp,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_ftp,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_ftp,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_ftp,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eononadcn_ftp,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eononadcn_ftp,'/D_EOnonADgtCN_pFWE*')); es1p=strcat(dir_grpcomp_eononadcn_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eononadcn_ftp,'/D_EOnonADltCN_pFWE*')); es2p=strcat(dir_grpcomp_eononadcn_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eononadcn_ftp,'/D_EOnonADgtCN_p0.1*')); es1p=strcat(dir_grpcomp_eononadcn_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eononadcn_ftp,'/D_EOnonADltCN_p0.1*')); es2p=strcat(dir_grpcomp_eononadcn_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eononadcn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% 4th: EOAD APOE+ vs. APOE- %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_grpcomp_eoadapoe=strcat(dir_grpcomp,'/EOAD_by_APOE_agesexmmse'); mkdir (dir_grpcomp_eoadapoe);
dir_grpcomp_eoadapoe_fbb=strcat(dir_grpcomp,'/EOAD_by_APOE_agesexmmse/fbb'); mkdir (dir_grpcomp_eoadapoe_fbb);
dir_grpcomp_eoadapoe_ftp=strcat(dir_grpcomp,'/EOAD_by_APOE_agesexmmse/ftp'); mkdir (dir_grpcomp_eoadapoe_ftp);

n_apoep=nnz(strcmp(dbm_eoad.APOEpos,'1'));
n_apoen=size(dbm_eoad,1)-n_apoep;

% Print a small log
            
Vars={'N APOE+'; 'N APOE-'; 'Covariates'};
Vals={n_apoep; n_apoen; 'Age,Sex,MMSE'};
tempT=table(Vars,Vals);
writetable(tempT, strcat(dir_grpcomp_eoadapoe,'/Log.csv'));

% FBB

apoep_wfbb=table2cell(dbm_eoad(1:n_apoep,{'Timepoint1_Path_FBB_wsuvr'}));
apoen_wfbb=table2cell(dbm_eoad(n_apoep+1:size(dbm_eoad,1),{'Timepoint1_Path_FBB_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eoadapoe_fbb);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = apoep_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = apoen_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eoad.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eoad.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c =  dbm_eoad.MMSE;
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'mmse';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'APOE+>APOE-';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'APOE+<APOE-';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'apoepMapoen_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'apoepMapoen_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'apoepLapoen_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'apoepLapoen_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_fbb,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_fbb,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_fbb,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_fbb,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eoadapoe_fbb,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eoadapoe_fbb,'/D_APOEgtAPOE-_pFWE*')); es1p=strcat(dir_grpcomp_eoadapoe_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eoadapoe_fbb,'/D_APOEgtAPOE-_pFWE*')); es2p=strcat(dir_grpcomp_eoadapoe_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eoadapoe_fbb,'/D_APOEgtAPOE-_p0.1*')); es1p=strcat(dir_grpcomp_eoadapoe_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eoadapoe_fbb,'/D_APOEgtAPOE-_p0.1*')); es2p=strcat(dir_grpcomp_eoadapoe_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eoadapoe)

% FTP

apoep_wftp=table2cell(dbm_eoad(1:n_apoep,{'Timepoint1_Path_FTP_wsuvr'}));
apoen_wftp=table2cell(dbm_eoad(n_apoep+1:size(dbm_eoad,1),{'Timepoint1_Path_FTP_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eoadapoe_ftp);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = apoep_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = apoen_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eoad.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eoad.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c =  dbm_eoad.MMSE;
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'mmse';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'APOE+>APOE-';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'APOE+<APOE-';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'apoepMapoen_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'apoepMapoen_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'apoepLapoen_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'apoepLapoen_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_ftp,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_ftp,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_ftp,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_ftp,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eoadapoe_ftp,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eoadapoe_ftp,'/D_APOEgtAPOE-_pFWE*')); es1p=strcat(dir_grpcomp_eoadapoe_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eoadapoe_ftp,'/D_APOEltAPOE-_pFWE*')); es2p=strcat(dir_grpcomp_eoadapoe_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eoadapoe_ftp,'/D_APOEgtAPOE-_p0.1*')); es1p=strcat(dir_grpcomp_eoadapoe_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eoadapoe_ftp,'/D_APOEltAPOE-_p0.1*')); es2p=strcat(dir_grpcomp_eoadapoe_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eoadapoe)


%%%%% No MMSE Version %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% 1st: EOAD vs. EOnonAD %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_grpcomp_eoadeononad=strcat(dir_grpcomp,'/EOADvsEOnonAD_agesexapoe'); mkdir (dir_grpcomp_eoadeononad);
dir_grpcomp_eoadeononad_fbb=strcat(dir_grpcomp,'/EOADvsEOnonAD_agesexapoe/fbb'); mkdir (dir_grpcomp_eoadeononad_fbb);
dir_grpcomp_eoadeononad_ftp=strcat(dir_grpcomp,'/EOADvsEOnonAD_agesexapoe/ftp'); mkdir (dir_grpcomp_eoadeononad_ftp);

n_eoad=nnz(strcmp(dbm_eoadeononad.CohortAssgn_x,'EOAD'));
n_eononad=size(dbm_eoadeononad,1)-n_eoad;

% Print a small log
            
Vars={'N EOAD'; 'N EOnonAD'; 'Covariates'};
Vals={n_eoad; n_eononad; 'Age,Sex,APOE'};
tempT=table(Vars,Vals);
writetable(tempT, strcat(dir_grpcomp_eoadeononad,'/Log.csv'));

% FBB

eoad_wfbb=table2cell(dbm_eoadeononad(1:n_eoad,{'Timepoint1_Path_FBB_wsuvr'}));
eononad_wfbb=table2cell(dbm_eoadeononad(n_eoad+1:size(dbm_eoadeononad,1),{'Timepoint1_Path_FBB_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eoadeononad_fbb);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = eoad_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = eononad_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eoadeononad.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eoadeononad.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c = str2double(dbm_eoadeononad.APOEpos);
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'apoe';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'EOAD>EOnonAD';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'EOAD<EOnonAD';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'eoadMeononad_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'eoadMeononad_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'eoadLeononad_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'eoadLeononad_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_fbb,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_fbb,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_fbb,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_fbb,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eoadeononad_fbb,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eoadeononad_fbb,'/D_EOADgtEOnonAD_pFWE*')); es1p=strcat(dir_grpcomp_eoadeononad_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eoadeononad_fbb,'/D_EOADltEOnonAD_pFWE*')); es2p=strcat(dir_grpcomp_eoadeononad_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eoadeononad_fbb,'/D_EOADgtEOnonAD_p0.1*')); es1p=strcat(dir_grpcomp_eoadeononad_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eoadeononad_fbb,'/D_EOADltEOnonAD_p0.1*')); es2p=strcat(dir_grpcomp_eoadeononad_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end


clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eoadeononad)

% FTP

eoad_wftp=table2cell(dbm_eoadeononad(1:n_eoad,{'Timepoint1_Path_FTP_wsuvr'}));
eononad_wftp=table2cell(dbm_eoadeononad(n_eoad+1:size(dbm_eoadeononad,1),{'Timepoint1_Path_FTP_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eoadeononad_ftp);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = eoad_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = eononad_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eoadeononad.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eoadeononad.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c = str2double(dbm_eoadeononad.APOEpos);
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'apoe';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'EOAD>EOnonAD';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'EOAD<EOnonAD';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'eoadMeononad_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'eoadMeononad_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'eoadLeononad_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'eoadLeononad_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_ftp,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_ftp,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_ftp,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadeononad_ftp,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eoadeononad_ftp,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eoadeononad_ftp,'/D_EOADgtEOnonAD_pFWE*')); es1p=strcat(dir_grpcomp_eoadeononad_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eoadeononad_ftp,'/D_EOADltEOnonAD_pFWE*')); es2p=strcat(dir_grpcomp_eoadeononad_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eoadeononad_ftp,'/D_EOADgtEOnonAD_p0.1*')); es1p=strcat(dir_grpcomp_eoadeononad_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eoadeononad_ftp,'/D_EOADltEOnonAD_p0.1*')); es2p=strcat(dir_grpcomp_eoadeononad_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadeononad_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% 2nd: EOAD vs. CN %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_grpcomp_eoadcn=strcat(dir_grpcomp,'/EOADvsCN_agesexapoe'); mkdir (dir_grpcomp_eoadcn);
dir_grpcomp_eoadcn_fbb=strcat(dir_grpcomp,'/EOADvsCN_agesexapoe/fbb'); mkdir (dir_grpcomp_eoadcn_fbb);
dir_grpcomp_eoadcn_ftp=strcat(dir_grpcomp,'/EOADvsCN_agesexapoe/ftp'); mkdir (dir_grpcomp_eoadcn_ftp);

n_eoad=nnz(strcmp(dbm_eoadcn.CohortAssgn_x,'EOAD'));
n_cn=size(dbm_eoadcn,1)-n_eoad;

% Print a small log
            
Vars={'N EOAD'; 'N CN'; 'Covariates'};
Vals={n_eoad; n_cn; 'Age,Sex,APOE'};
tempT=table(Vars,Vals);
writetable(tempT, strcat(dir_grpcomp_eoadcn,'/Log.csv'));

% FBB

eoad_wfbb=table2cell(dbm_eoadcn(1:n_eoad,{'Timepoint1_Path_FBB_wsuvr'}));
cn_wfbb=table2cell(dbm_eoadcn(n_eoad+1:size(dbm_eoadcn,1),{'Timepoint1_Path_FBB_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eoadcn_fbb);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = eoad_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = cn_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eoadcn.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eoadcn.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c = str2double(dbm_eoadcn.APOEpos);
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'apoe';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'EOAD>CN';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'EOAD<CN';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'eoadMcn_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'eoadMcn_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'eoadLcn_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'eoadLcn_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_fbb,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_fbb,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_fbb,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_fbb,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eoadcn_fbb,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eoadcn_fbb,'/D_EOADgtCN_pFWE*')); es1p=strcat(dir_grpcomp_eoadcn_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eoadcn_fbb,'/D_EOADltCN_pFWE*')); es2p=strcat(dir_grpcomp_eoadcn_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eoadcn_fbb,'/D_EOADgtCN_p0.1*')); es1p=strcat(dir_grpcomp_eoadcn_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eoadcn_fbb,'/D_EOADltCN_p0.1*')); es2p=strcat(dir_grpcomp_eoadcn_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eoadcn)

% FTP

eoad_wftp=table2cell(dbm_eoadcn(1:n_eoad,{'Timepoint1_Path_FTP_wsuvr'}));
cn_wftp=table2cell(dbm_eoadcn(n_eoad+1:size(dbm_eoadcn,1),{'Timepoint1_Path_FTP_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eoadcn_ftp);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = eoad_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = cn_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eoadcn.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eoadcn.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c = str2double(dbm_eoadcn.APOEpos);
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'apoe';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'EOAD>CN';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'EOAD<CN';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'eoadMcn_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'eoadMcn_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'eoadLcn_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'eoadLcn_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_ftp,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_ftp,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_ftp,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadcn_ftp,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eoadcn_ftp,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eoadcn_ftp,'/D_EOADgtCN_pFWE*')); es1p=strcat(dir_grpcomp_eoadcn_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eoadcn_ftp,'/D_EOADltCN_pFWE*')); es2p=strcat(dir_grpcomp_eoadcn_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eoadcn_ftp,'/D_EOADgtCN_p0.1*')); es1p=strcat(dir_grpcomp_eoadcn_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eoadcn_ftp,'/D_EOADltCN_p0.1*')); es2p=strcat(dir_grpcomp_eoadcn_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eoadcn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% 3rd: EOnonAD vs. CN %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_grpcomp_eononadcn=strcat(dir_grpcomp,'/EOnonADvsCN_agesexapoe'); mkdir (dir_grpcomp_eononadcn);
dir_grpcomp_eononadcn_fbb=strcat(dir_grpcomp,'/EOnonADvsCN_agesexapoe/fbb'); mkdir (dir_grpcomp_eononadcn_fbb);
dir_grpcomp_eononadcn_ftp=strcat(dir_grpcomp,'/EOnonADvsCN_agesexapoe/ftp'); mkdir (dir_grpcomp_eononadcn_ftp);

n_eononad=nnz(strcmp(dbm_eononadcn.CohortAssgn_x,'EOnonAD'));
n_cn=size(dbm_eononadcn,1)-n_eononad;

% Print a small log
            
Vars={'N EOnonAD'; 'N CN'; 'Covariates'};
Vals={n_eononad; n_cn; 'Age,Sex,APOE'};
tempT=table(Vars,Vals);
writetable(tempT, strcat(dir_grpcomp_eononadcn,'/Log.csv'));

% FBB

eononad_wfbb=table2cell(dbm_eononadcn(1:n_eononad,{'Timepoint1_Path_FBB_wsuvr'}));
cn_wfbb=table2cell(dbm_eononadcn(n_eononad+1:size(dbm_eononadcn,1),{'Timepoint1_Path_FBB_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eononadcn_fbb);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = eononad_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = cn_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eononadcn.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eononadcn.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c = str2double(dbm_eononadcn.APOEpos);
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'apoe';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'EOnonAD>CN';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'EOnonAD<CN';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'eononadMcn_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'eononadMcn_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'eononadLcn_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'eononadLcn_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_fbb,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_fbb,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_fbb,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_fbb,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eononadcn_fbb,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eononadcn_fbb,'/D_EOnonADgtCN_pFWE*')); es1p=strcat(dir_grpcomp_eononadcn_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eononadcn_fbb,'/D_EOnonADltCN_pFWE*')); es2p=strcat(dir_grpcomp_eononadcn_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eononadcn_fbb,'/D_EOnonADgtCN_p0.1*')); es1p=strcat(dir_grpcomp_eononadcn_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eononadcn_fbb,'/D_EOnonADltCN_p0.1*')); es2p=strcat(dir_grpcomp_eononadcn_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eononadcn)

% FTP

eononad_wftp=table2cell(dbm_eononadcn(1:n_eononad,{'Timepoint1_Path_FTP_wsuvr'}));
cn_wftp=table2cell(dbm_eononadcn(n_eononad+1:size(dbm_eononadcn,1),{'Timepoint1_Path_FTP_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eononadcn_ftp);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = eononad_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = cn_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eononadcn.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eononadcn.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.cov(3).c = str2double(dbm_eononadcn.APOEpos);
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'apoe';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'EOnonAD>CN';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'EOnonAD<CN';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'eononadMcn_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'eononadMcn_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'eononadLcn_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'eononadLcn_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_ftp,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_ftp,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_ftp,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eononadcn_ftp,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eononadcn_ftp,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eononadcn_ftp,'/D_EOnonADgtCN_pFWE*')); es1p=strcat(dir_grpcomp_eononadcn_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eononadcn_ftp,'/D_EOnonADltCN_pFWE*')); es2p=strcat(dir_grpcomp_eononadcn_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eononadcn_ftp,'/D_EOnonADgtCN_p0.1*')); es1p=strcat(dir_grpcomp_eononadcn_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eononadcn_ftp,'/D_EOnonADltCN_p0.1*')); es2p=strcat(dir_grpcomp_eononadcn_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eononadcn_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eononadcn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% 4th: EOAD APOE+ vs. APOE- %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_grpcomp_eoadapoe=strcat(dir_grpcomp,'/EOAD_by_APOE_agesex'); mkdir (dir_grpcomp_eoadapoe);
dir_grpcomp_eoadapoe_fbb=strcat(dir_grpcomp,'/EOAD_by_APOE_agesex/fbb'); mkdir (dir_grpcomp_eoadapoe_fbb);
dir_grpcomp_eoadapoe_ftp=strcat(dir_grpcomp,'/EOAD_by_APOE_agesex/ftp'); mkdir (dir_grpcomp_eoadapoe_ftp);

n_apoep=nnz(strcmp(dbm_eoad.APOEpos,'1'));
n_apoen=size(dbm_eoad,1)-n_apoep;

% Print a small log
            
Vars={'N APOE+'; 'N APOE-'; 'Covariates'};
Vals={n_apoep; n_apoen; 'Age,Sex'};
tempT=table(Vars,Vals);
writetable(tempT, strcat(dir_grpcomp_eoadapoe,'/Log.csv'));

% FBB

apoep_wfbb=table2cell(dbm_eoad(1:n_apoep,{'Timepoint1_Path_FBB_wsuvr'}));
apoen_wfbb=table2cell(dbm_eoad(n_apoep+1:size(dbm_eoad,1),{'Timepoint1_Path_FBB_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eoadapoe_fbb);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = apoep_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = apoen_wfbb;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eoad.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eoad.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'APOE+>APOE-';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'APOE+<APOE-';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'apoepMapoen_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'apoepMapoen_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'apoepLapoen_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'apoepLapoen_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_fbb,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_fbb,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_fbb,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_fbb,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eoadapoe_fbb,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eoadapoe_fbb,'/D_APOEgtAPOE-_pFWE*')); es1p=strcat(dir_grpcomp_eoadapoe_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eoadapoe_fbb,'/D_APOEgtAPOE-_pFWE*')); es2p=strcat(dir_grpcomp_eoadapoe_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eoadapoe_fbb,'/D_APOEgtAPOE-_p0.1*')); es1p=strcat(dir_grpcomp_eoadapoe_fbb,'/',es1.name);

if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
end

es2=dir(strcat(dir_grpcomp_eoadapoe_fbb,'/D_APOEgtAPOE-_p0.1*')); es2p=strcat(dir_grpcomp_eoadapoe_fbb,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_fbb,'/3DRend_FBB_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FBB_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eoadapoe)

% FTP

apoep_wftp=table2cell(dbm_eoad(1:n_apoep,{'Timepoint1_Path_FTP_wsuvr'}));
apoen_wftp=table2cell(dbm_eoad(n_apoep+1:size(dbm_eoad,1),{'Timepoint1_Path_FTP_wsuvr'}));

spm('defaults','PET');
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(dir_grpcomp_eoadapoe_ftp);
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = apoep_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = apoen_wftp;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c =  dbm_eoad.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c =  double(dbm_eoad.sex_rcd);
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = cellstr(emask);
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'APOE+>APOE-';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'APOE+<APOE-';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.tspm.basename = 'apoepMapoen_thresfwe05';
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.titlestr = '';
matlabbatch{5}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{5}.spm.stats.results.conspec.extent = 100;
matlabbatch{5}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.export{1}.tspm.basename = 'apoepMapoen_thresunc001';
matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.titlestr = '';
matlabbatch{6}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{6}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{6}.spm.stats.results.units = 1;
matlabbatch{6}.spm.stats.results.export{1}.tspm.basename = 'apoepLapoen_thresfwe05';
matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
matlabbatch{7}.spm.stats.results.conspec.contrasts = 2;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{7}.spm.stats.results.conspec.extent = 100;
matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{7}.spm.stats.results.units = 1;
matlabbatch{7}.spm.stats.results.export{1}.tspm.basename = 'apoepLapoen_thresunc001';
matlabbatch{8}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_ftp,'/spmT_0001.nii'));
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{8}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{8}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{9}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_ftp,'/spmT_0001.nii'));
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{9}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{9}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{10}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_ftp,'/spmT_0002.nii'));
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.threshdesc.fwe.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{10}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{10}.spm.tools.cat.tools.T2x.atlas = 'None';
matlabbatch{11}.spm.tools.cat.tools.T2x.data_T2x = cellstr(strcat(dir_grpcomp_eoadapoe_ftp,'/spmT_0002.nii'));
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.sel = 4;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.threshdesc.uncorr.thresh001 = 0.001;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.inverse = 0;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.thresh05 = 0.05;
matlabbatch{11}.spm.tools.cat.tools.T2x.conversion.cluster.fwe2.noniso = 1;
matlabbatch{11}.spm.tools.cat.tools.T2x.atlas = 'None';

batch_str=strcat(dir_grpcomp_eoadapoe_ftp,'/Batch');
            save(char(batch_str),'matlabbatch'); %saves batch for user
			spm_jobman('run',matlabbatch);

% Get Automated 3D renderings 

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
es1=dir(strcat(dir_grpcomp_eoadapoe_ftp,'/D_APOEgtAPOE-_pFWE*')); es1p=strcat(dir_grpcomp_eoadapoe_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eoadapoe_ftp,'/D_APOEltAPOE-_pFWE*')); es2p=strcat(dir_grpcomp_eoadapoe_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear es1 es2

es1=dir(strcat(dir_grpcomp_eoadapoe_ftp,'/D_APOEgtAPOE-_p0.1*')); es1p=strcat(dir_grpcomp_eoadapoe_ftp,'/',es1.name);
if size(es1,1)==1
rend_img=es1p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

es2=dir(strcat(dir_grpcomp_eoadapoe_ftp,'/D_APOEltAPOE-_p0.1*')); es2p=strcat(dir_grpcomp_eoadapoe_ftp,'/',es2.name);
if size(es2,1)==1
rend_img=es2p;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(dir_grpcomp_eoadapoe_ftp,'/3DRend_FTP_',ff,'.jpg');
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Options_BrainNet_3D_FTP_EffectSize_unc.mat';
BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
end

clear matlabbatch Vars Vals tempT es1 es2
cd (dir_grpcomp_eoadapoe)

%% end

cd (path_master)
clear;

