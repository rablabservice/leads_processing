function outfiles = save_mask_infcblgm(aparcf, suitf, in_dir, out_dir, overwrite)
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
    % suitf : char or str array
    %   Path to the SUIT atlas in native MRI space
    % in_dir : char or str array
    %   The input directory. If aparcf is empty, this is where the
    %   function looks for the aparc_dat+aseg.nii file. This parameter is
    %   disregarded if aparcf is not empty
    % out_dir : char or str array
    %   The output directory. If out_dir is empty, the mask is saved in
    %   the same directory as the aparc_dat+aseg.nii file
    % overwrite : logical, optional
    %   If true, overwrite existing file
    %
    % Files created
    % -------------
    % - <scan_tag>_cbl-suit.nii
    % - <scan_tag>_mask-cblgm.nii
    % - <scan_tag>_mask-infcblgm.nii
    % ------------------------------------------------------------------
    arguments
        aparcf {mustBeText} = ''
        suitf {mustBeText} = ''
        in_dir {mustBeText} = ''
        out_dir {mustBeText} = ''
        overwrite logical = false
    end

    function smooth_mask(infile)
        % Smooth infile by 8mm^3 FWHM and save as s<infile>
        prefix = 's';
        clear matlabbatch;
        matlabbatch{1}.spm.spatial.smooth.data = cellstr(infile);
        matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
        matlabbatch{1}.spm.spatial.smooth.dtype = spm_type('float32');
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = prefix;
        spm_jobman('run', matlabbatch);
        fprintf('  * Saved %s\n', basename(add_presuf(infile, prefix)));
    end

    % Define aparc_dat indices for cerebellar gray matter
    mask_idx = [8; 47];

    % Format inputs
    [aparcf, out_dir] = format_mask_inputs(aparcf, in_dir, out_dir);
    scan_tag = get_scan_tag(aparcf);
    mustBeFile(aparcf);
    mustBeFile(suitf);

    % Define outputs
    outfiles.mask_cblgm = fullfile(out_dir, append(scan_tag, '_mask-cblgm.nii'));
    % outfiles.cbl_suit = fullfile(out_dir, append(scan_tag, '_cbl-suit.nii'));
    outfiles.mask_infcblgm = fullfile(out_dir, append(scan_tag, '_mask-infcblgm.nii'));

    % Save the mask
    mask_cblgm = nii_labels_to_mask(aparcf, mask_idx, outfiles.mask_cblgm, overwrite);

    % % Inverse warp the SUIT template to native MRI space
    % if isfile(outfiles.cbl_suit) && ~overwrite
    %     fprintf('  * Reverse-normalized SUIT template already exists, will not overwrite\n');
    % else
    %     template_suitf = '/mnt/coredata/Projects/Resources/scripts/rCerebellum-SUIT.nii';
    %     mustBeFile(template_suitf);
    %     interp = 0;
    %     vox = [1 1 1];
    %     prefix = 'v';
    %     bb = [Inf Inf Inf; Inf Inf Inf];
    %     outf = apply_invwarp_from_mni(template_suitf, iyf, interp, vox, prefix, bb, overwrite);
    %     movefile(outf, outfiles.cbl_suit);
    %     fprintf('  * Saved %s\n', basename(outfiles.cbl_suit));
    % end

    % Define intermediary files
    interfs.mask_keep = fullfile(out_dir, append(scan_tag, '_mask-keep.nii'));
    interfs.smask_keep = add_presuf(interfs.mask_keep, 's');
    interfs.mask_toss = fullfile(out_dir, append(scan_tag, '_mask-toss.nii'));
    interfs.smask_toss = add_presuf(interfs.mask_toss, 's');
    interfs.keep_gt_toss = fullfile(out_dir, 'keep_gt_toss.nii');

    % Create a mask of SUIT subregions we want to keep, then smooth it
    labels_keep = [6, 8:28, 33:34]';
    mask_keep = nii_labels_to_mask(suitf, labels_keep, interfs.mask_keep, overwrite);
    smooth_mask(interfs.mask_keep);
    smask_keep = spm_read_vols(spm_vol(interfs.smask_keep));

    % Create a mask of SUIT subregions we want to keep, then smooth it
    labels_toss = [0:5, 7]';
    mask_toss = nii_labels_to_mask(suitf, labels_toss, interfs.mask_toss, overwrite);
    smooth_mask(interfs.mask_toss);
    smask_toss = spm_read_vols(spm_vol(interfs.smask_toss));

    % Save the inferior cerebellar gray matter mask
    keep_gt_toss = smask_keep > smask_toss;
    out_img = spm_vol(suitf);
    out_img.fname = interfs.keep_gt_toss;
    out_img.dt = [spm_type('uint8') 0];
    spm_write_vol(out_img, keep_gt_toss);

    mask_infcblgm = mask_cblgm & (smask_keep > smask_toss);
    out_img = spm_vol(aparcf);
    out_img.fname = outfiles.mask_infcblgm;
    out_img.dt = [spm_type('uint8') 0];
    spm_write_vol(out_img, mask_infcblgm);
    fprintf('  * Saved %s\n', basename(outfiles.mask_infcblgm));

    % Delete intermediary files
    % structfun(@delete, interfs);
end
