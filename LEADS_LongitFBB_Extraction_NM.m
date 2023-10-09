% Temporary script to calculate longitudinal FBB SUVRs 
% Calculating using a straight average between whole cerebellum, brainstem,
% and eroded wm mask 
% June 2023; nidhi.mundada@ucsf.edu
% Used primarily for Ganna's Amyloid Chronicity Project

path_processed='/mnt/coredata/Projects/LEADS/data_f7p1/processed/';

listfbbs=dir(strcat(path_processed,'LDS*/*/FBB*/rLDS*FBB*.nii'));

% create empty matrix in which to store values
M_bigref_sf = zeros(size(listfbbs,1),1); 
M_globalscore=zeros(size(listfbbs,1),1);
M_sz=zeros(size(listfbbs,1),6);

for i=1:size(listfbbs,1)
    fbbdir=listfbbs(i,1).folder;
    fbbfile=strcat(fbbdir,'/',listfbbs(i,1).name);
    
    % read rFBB 
    img=spm_vol(fbbfile); % Reading img header from the loop
    img1=spm_read_vols(img); % Reading img values 
    
    % read whole cerebellum and extract mean value from rFBB*.nii 
    wholecere=spm_vol(char(strcat(fbbdir,'/wholecbl_ref_mask.nii')));
    wholecere1=spm_read_vols(wholecere);

    img_mask_wholecere=img1.*wholecere1; % Creating the masked image
    img_mask_wholecere=nonzeros(img_mask_wholecere);
    img_mask_wholecere=img_mask_wholecere(~isnan(img_mask_wholecere));
    wmean_wholecere=mean(img_mask_wholecere);

    % read eroded wm and extract mean value from rFBB*.nii 
    erodedwm=dir(strcat(fbbdir,'/s*raparc+aseg_wm_thr0p7.nii')); 
    erodedwm=spm_vol(char(strcat(erodedwm.folder,'/', erodedwm.name)));
    erodedwm1=spm_read_vols(erodedwm);

    img_mask_erodedwm=img1.*erodedwm1; % Creating the masked image
    img_mask_erodedwm=nonzeros(img_mask_erodedwm);
    img_mask_erodedwm=img_mask_erodedwm(~isnan(img_mask_erodedwm));
    wmean_erodedwm=mean(img_mask_erodedwm);
    
    % read bs and extract mean value from rFBB*.nii 
    bs=dir(strcat(fbbdir,'/LDS*raparc+aseg_bs.nii')); 
    bs=spm_vol(char(strcat(bs.folder,'/', bs.name)));
    bs1=spm_read_vols(bs);

    img_mask_bs=img1.*bs1; % Creating the masked image
    img_mask_bs=nonzeros(img_mask_bs);
    img_mask_bs=img_mask_bs(~isnan(img_mask_bs));
    wmean_bs=mean(img_mask_bs);
    
    bigref=[wmean_wholecere,wmean_erodedwm,wmean_bs];
    bigref_sf=mean(bigref);
    M_bigref_sf(i,1)=bigref_sf;
    
    % extract global score from rFBB*.nii
    aparcscan=dir(strcat(fbbdir,'/LDS*raparc+aseg.nii'));
    aparcscan=strcat(aparcscan.folder,'/',aparcscan.name);

    rfbb=spm_vol(fbbfile);
    Vaparc=spm_vol(aparcscan);

    rfbb_vol=spm_read_vols(rfbb);
    [sz1,sz2,sz3]=size(rfbb_vol);
    rrfbb_vol=reshape(rfbb_vol,sz1*sz2*sz3,1);

    aparc=spm_read_vols(Vaparc);
    rparc=reshape(aparc,sz1*sz2*sz3,1);

    ind2=find((rparc==1003 | rparc==1012 | rparc==1014 ...
        | rparc==1018 | rparc==1019 | rparc==1020 ...
        | rparc==1027 | rparc==1028 | rparc==1032 ...
        | rparc==1009 | rparc==1015 | rparc==1030 ...
        | rparc==1008 | rparc==1025 | rparc==1029 ...
        | rparc==1031 | rparc==1002 | rparc==1023 | rparc==1010 ...
        | rparc==1026 | rparc==2003 | rparc==2012 | rparc==2014 ...
        | rparc==2018 | rparc==2019 | rparc==2020 ...
        | rparc==2027 | rparc==2028 | rparc==2032 ...
        | rparc==2009 | rparc==2015 | rparc==2030 ...
        | rparc==2008 | rparc==2025 | rparc==2029 ...
        | rparc==2031 | rparc==2002 | rparc==2023 | rparc==2010 ...
        | rparc==2026) & rrfbb_vol>0);

    compsuvr_t2=mean(rrfbb_vol(ind2));
    M_globalscore(i,1)=compsuvr_t2;
    
    %size 
    M_sz(i,1)=wmean_wholecere;
    M_sz(i,2)=wmean_erodedwm;
    M_sz(i,3)=wmean_bs;
    M_sz(i,4)=size(img_mask_wholecere,1);
    M_sz(i,5)=size(img_mask_erodedwm,1);
    M_sz(i,6)=size(img_mask_bs,1);
    
    clearvars -except path_processed listfbbs M_bigref_sf M_globalscore M_sz
    
end 

T_bigref_sf=array2table(M_bigref_sf);
T_globalscore=array2table(M_globalscore);
T_bigref_sf.Properties.VariableNames={'ScalingFactor_CompWM'};
T_globalscore.Properties.VariableNames={'rFBB_GlobalScore'};

T_sz=array2table(M_sz);
T_sz.Properties.VariableNames={'Wholecere_Step4', 'Erodedwm_Step4','BS_Step4','Wholecere_Step4_Size', 'Erodedwm_Step4_Size','BS_Step4_Size'};

x=struct2table(listfbbs);
x=x(:,1);
x.Properties.VariableNames={'Filename'};

db=[x T_bigref_sf T_globalscore T_sz];

filename = sprintf('/mnt/coredata/Projects/LEADS/data_f7p1/extraction/FBB_Modified_CompWM_Extraction_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
writetable(db,filename,'WriteRowNames',true)

