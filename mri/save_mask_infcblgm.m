function outfiles = save_mask_infcblgm(aparcf, iyf, in_dir, out_dir, overwrite, verbose)
    % Save the inferior cerebellar gray matter mask
    %
    % Usage
    % -----
    % >> save_mask_infcblgm(aparcf)
    %
    % Parameters
    % ----------
    % aparcf : char or str array
    %   Path to the aparc+aseg.nii file
    % iyf : char or str array
    %   Path to the inverse warp file
    % in_dir : char or str array
    %   The input directory. If aparcf is empty, this is where the
    %   function looks for the aparc+aseg.nii file. This parameter is
    %   disregarded if aparcf is not empty
    % out_dir : char or str array
    %   The output directory. If out_dir is empty, the mask is saved in
    %   the same directory as the aparc+aseg.nii file
    % overwrite : logical, optional
    %   If true, overwrite existing file
    % verbose : logical, optional
    %   If true, print diagnostic information
    %
    % Files created
    % -------------
    % - <out_dir>/<scan_tag>_cbl-suit.nii
    % - <out_dir>/<scan_tag>_mask-infcblgm.nii
    % ------------------------------------------------------------------
    arguments
        aparcf {mustBeText} = ''
        in_dir {mustBeText} = ''
        out_dir {mustBeText} = ''
        overwrite logical = false
        verbose logical = true
    end

    % Define aparc indices for cerebellar gray matter
    mask_idx = [8; 47];

    % Format inputs
    [aparcf, out_dir] = format_mask_inputs(aparcf, in_dir, out_dir);
    scan_tag = get_scan_tag(aparcf);

    % Save the mask
    cblgmf = fullfile(out_dir, append(scan_tag, '_mask-cblgm.nii'));
    nii_labels_to_mask(aparcf, mask_idx, outfile, overwrite, verbose);

    % Inverse warp the SUIT template to subject MRI space
    template_suitf = '/mnt/coredata/Projects/Resources/scripts/rCerebellum-SUIT.nii'

    spm('defaults','PET');
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(iyf);
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(template_suitf);
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

    % Code straight from ADNI processing sent from Deniz and available at /home/jagust/dkorman/matlab/MakeSUITCerebMask_ADNI.m

    Vaparc=spm_vol(aparcf);
    aparc=spm_read_vols(Vaparc);
    [sz1,sz2,sz3]=size(aparc);
    raparc=reshape(aparc,sz1*sz2*sz3,1);

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
    ind=find((raparc==8 | raparc==47) & rskeep>rstoss);
    szwholecere=length(ind);
    rcereaparc=zeros(sz1*sz2*sz3,1);
    rcereaparc(ind)=ones(length(ind),1);

    cereaparc=reshape(rcereaparc,sz1,sz2,sz3);
    Vcereaparc=Vaparc;
    Vcereaparc.fname=[pathftp '/infcblg_ref_mask.nii']; %% saved the inferior cbl gray mask
    spm_write_vol(Vcereaparc,cereaparc);

    % Create the inferior cerebellar gray matter mask


end
