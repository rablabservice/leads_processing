%%% script to run when in the FBB folder where an FBB needs to be processed
%%% outputs: Global FBB, Mask ref region whole CBL, SUVR image
%%% inputs FBB scan, nu MRI, aparc+aseg MRI

%nii_setOrigin(fbbscan); % reset origin
% Different module to set origin
file = deblank(fbbscan);
st.vol = spm_vol(file);
vs = st.vol.mat\eye(4);
vs(1:3,4) = (st.vol.dim+1)/2;
spm_get_space(st.vol.fname,inv(vs));

% coreg scan to respective mri
spm('defaults','PET');
clear matlabbatch;
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(nuscan);
matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(fbbscan);
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
clear matlabbatch;

r1n=spm_vol(aparcscan); % reading the aparc-aseg
r1=spm_read_vols(r1n);

img=spm_vol(rfbbscan); % Reading the coregistered fbb scan we just estimated
img1=spm_read_vols(img);

%%% Extract ref region %%%
ind_cbl=[7 8 46 47]';
Mcbl = zeros(1,size(ind_cbl,1)); % create empty matrix in which to store values from all the ROIs
Mcbl_sz = zeros(1,size(ind_cbl,1));
for f = 1:size(ind_cbl,1) % Looping through the ROIs
    reg = ind_cbl(f);  % Selecting the ROI from the loop
    mask = (r1 == reg); % Creating the logical mask to separate the ROI from the loop
    img_mask=img1.*mask; % Creating the masked image
    temp=nonzeros(img_mask);
    temp=temp(~isnan(temp));
    ext_val=mean(temp);
    vec(f)=ext_val; % Store values for the individual ROIs across the images
    vec_sz(f)=size(temp,1);
end
Mcbl(1,:)=vec; % save the values ROI-wise (column-wise)
Mcbl_sz(1,:)=vec_sz;  % save the roi size
wmean_cbl = sum(Mcbl_sz.*Mcbl,2)./sum(Mcbl_sz,2);

%% let's save the mask for the reference region and store the SUVR image
exp=char(strcat('i1/',num2str(wmean_cbl)));
newfname=char(strcat(fbbscan(1:end-4),'_suvr_cbl.nii'));

spm('defaults','PET');
clear matlabbatch;
matlabbatch{1}.spm.util.imcalc.input = cellstr(aparcscan);
matlabbatch{1}.spm.util.imcalc.output = 'wholecbl_ref_mask';
matlabbatch{1}.spm.util.imcalc.outdir = cellstr(pathfbb);
matlabbatch{1}.spm.util.imcalc.expression = '(i1==7) | ( i1==8) | (i1==46) | (i1==47)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 0;
matlabbatch{1}.spm.util.imcalc.options.dtype = spm_type('uint8');
matlabbatch{2}.spm.util.imcalc.input = cellstr(rfbbscan);
matlabbatch{2}.spm.util.imcalc.output = newfname;
matlabbatch{2}.spm.util.imcalc.outdir = {''};
matlabbatch{2}.spm.util.imcalc.expression = exp;
matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{2}.spm.util.imcalc.options.mask = 0;
matlabbatch{2}.spm.util.imcalc.options.interp = 1;
matlabbatch{2}.spm.util.imcalc.options.dtype = spm_type('float32');
spm_jobman('run',matlabbatch);
clear matlabbatch;

% New module to extract a global score with the updated ADNI method
% (weighted average across FS7.1 regions)
Vsuvr=spm_vol(newfname);
Vaparc=spm_vol(aparcscan);

suvr=spm_read_vols(Vsuvr);
[sz1,sz2,sz3]=size(suvr);
rsuvr=reshape(suvr,sz1*sz2*sz3,1);

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
    | rparc==2026) & rsuvr>0);

compsuvr_t2=mean(rsuvr(ind2));
csvwrite(strcat(pathfbb,'/adni-global-amyloid-index_suvr-wcbl_fs7-1.csv'),compsuvr_t2); %%  Saved the global summary score as a csv in the FBB folder
