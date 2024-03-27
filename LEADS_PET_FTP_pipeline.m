%%% script to run when in the FTP folder where an FTP needs to be processed
%%% outputs: Mask ref region inferior CBL gray, SUVR image
%%% inputs FTP scan, nu MRI, aparc+aseg MRI, SUIT Template

%nii_setOrigin(ftpscan);
% Different module to set origin
file = deblank(ftpscan);
st.vol = spm_vol(file);
vs = st.vol.mat\eye(4);
vs(1:3,4) = (st.vol.dim+1)/2;
spm_get_space(st.vol.fname,inv(vs));

spm('defaults','PET');
clear matlabbatch;
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(nuscan);
matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(ftpscan);
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch); clear matlabbatch;

%% Working to get the SUIT template reverse normalized in the native subject space

copyfile('/mnt/coredata/Projects/LEADS/script_f7p1/templates/rCerebellum-SUIT.nii',pathftp)

spm('defaults','PET');
clear matlabbatch;
matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(strcat(pathftp,'/y_revnorm.nii'));
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(strcat(pathftp,'/rCerebellum-SUIT.nii'));
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
matlabbatch{2}.spm.spatial.coreg.write.ref = cellstr(nuscan);
matlabbatch{2}.spm.spatial.coreg.write.source = cellstr(strcat(pathftp,'/wrCerebellum-SUIT.nii'));
matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r';
spm_jobman('run',matlabbatch); clear matlabbatch;

%%% we have the reverse normalized SUIT template now.
%%% Let's take into account the aparc+aseg

%%% Code straight from ADNI processing sent from Deniz and available at /home/jagust/dkorman/matlab/MakeSUITCerebMask_ADNI.m %%%

Vaparc=spm_vol(aparcscan);
aparc=spm_read_vols(Vaparc);
[sz1,sz2,sz3]=size(aparc);
raparc=reshape(aparc,sz1*sz2*sz3,1);

Vsuvr=spm_vol(rftpscan);
suvr=spm_read_vols(Vsuvr);
rsuvr=reshape(suvr,sz1*sz2*sz3,1);

Vcere=spm_vol(char(strcat(pathftp,'/rwrCerebellum-SUIT.nii')));
cere=spm_read_vols(Vcere);
rcere=reshape(cere,sz1*sz2*sz3,1);
% find voxels we want to keep and toss in reverse-normalized cerebellum
% atlas
indkeep=find(rcere==6 | (rcere>=8 & rcere<=28) | rcere==33 | rcere==34);
indtoss=find(rcere<=5 | rcere==7);
rkeep=zeros(sz1*sz2*sz3,1);
rtoss=zeros(sz1*sz2*sz3,1);
% create binary masks for voxels we want to keep or toss
rkeep(indkeep)=ones(length(indkeep),1);
rtoss(indtoss)=ones(length(indtoss),1);
keep=reshape(rkeep,sz1,sz2,sz3);
toss=reshape(rtoss,sz1,sz2,sz3);
skeep=zeros(sz1,sz2,sz3);
stoss=zeros(sz1,sz2,sz3);
% smooth the binary masks for voxels we want to keep or toss, doing this bc
% there is not perfect overlap bw freesurfer's gray matter segmentation of
% cerebellum and the reversenormalized mask, so want a freesurfer gray
% matter voxel to be characterized in keep or toss group depending on how
% close it is to keep or toss regions in reverse normalized cerebellum
% template, even if it isn't defined as part of cerebellum in the template
spm_smooth(keep,skeep,[8 8 8]);
spm_smooth(toss,stoss,[8 8 8]);
rskeep=reshape(skeep,sz1*sz2*sz3,1);
rstoss=reshape(stoss,sz1*sz2*sz3,1);
ind=find((raparc==8 | raparc==47) & rsuvr>0 & rskeep>rstoss);
szwholecere=length(ind);
rcereaparc=zeros(sz1*sz2*sz3,1);
rcereaparc(ind)=ones(length(ind),1);

cereaparc=reshape(rcereaparc,sz1,sz2,sz3);
Vcereaparc=Vaparc;
Vcereaparc.fname=[pathftp '/infcblg_ref_mask.nii']; %% saved the inferior cbl gray mask
spm_write_vol(Vcereaparc,cereaparc);

meaninfcere=mean(rsuvr(ind)); %% stored the ref region value to create the SUVR image

exp=char(strcat('i1/',num2str(meaninfcere)));
newfname=char(strcat(ftpscan(1:end-4),'_suvr_infcblg.nii'));

spm('defaults','PET');
clear matlabbatch;
matlabbatch{1}.spm.util.imcalc.input = cellstr(rftpscan);
matlabbatch{1}.spm.util.imcalc.output = newfname;
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = exp;
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = spm_type('float32');
spm_jobman('run',matlabbatch); clear matlabbatch;
