%% Service to create 3D Renderings for Wmap %%
%%% Sep 2021: Slight change

%%%%% code to create 3D surface renderings %%%%

rend_bg='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
rend_img=tempwmapname_gm;
[pp,ff,ee]=spm_fileparts(rend_img);
rend_fname=strcat(pp,'/3DRend_',ff,'.jpg');

if strcmp(tempmod,'FTP')==1
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Option_BrainNet_Wmap_FTP.mat';
elseif strcmp(tempmod,'FBB')==1
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Option_BrainNet_Wmap_FBB.mat';
elseif strcmp(tempmod,'FDG')==1 % FDG and MRI have the same bounds i.e. -1.65 and -5
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Option_BrainNet_Wmap_MRI.mat';
elseif strcmp(tempmod,'MRI')==1
rend_opt='/mnt/coredata/Projects/LEADS/script_f7p1/service/BrainNetViewer/LEADS_Options/Option_BrainNet_Wmap_MRI.mat';
end

BrainNet_MapCfg(rend_bg,rend_img,rend_opt,rend_fname);
close all
