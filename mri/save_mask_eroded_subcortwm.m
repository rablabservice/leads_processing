function outfiles = save_mask_eroded_subcortwm( ...
    aparcf, fid, overwrite, smooth_by, erosion_thresh, in_dir, out_dir ...
)
    % Load the aparc+aseg and save eroded subcortical white matter mask
    %
    % Parameters
    % ----------
    % aparcf : char or str array
    %   Path to the aparc+aseg.nii file
    % fid : int, optional
    %   File identifier for logging (default is 1 for stdout)
    % overwrite : logical, optional
    %   If true, overwrite existing file
    % smooth_by : double
    %   FWHM of the smoothing kernel applied to the subcortical WM mask.
    %   Default is 8mm.
    % erosion_thresh : double
    %   Threshold applied to the smoothed mask to get the eroded mask.
    %   Values > erosion_thresh are retained. Default is 0.7.
    % in_dir : char or str array
    %   The input directory. If aparcf is empty, this is where the
    %   function looks for the aparc+aseg.nii file. This parameter is
    %   disregarded if aparcf is not empty
    % out_dir : char or str array
    %   The output directory. If out_dir is empty, the mask is saved in
    %   the same directory as the aparc+aseg.nii file
    %
    % Files created
    % -------------
    % - <scan_tag>_mask-eroded-subcortwm.nii
    % - <scan_tag>_mask-subcortwm.nii
    % - s<scan_tag>_mask-subcortwm.nii
    % ------------------------------------------------------------------
    arguments
        aparcf {mustBeText} = ''
        fid {mustBeNumeric} = 1
        overwrite logical = false
        smooth_by {mustBeNumeric} = 8
        erosion_thresh {mustBeNumeric} = 0.7
        in_dir {mustBeText} = ''
        out_dir {mustBeText} = ''
    end

    % Define aparc indices for the subcortical white matter
    mask_idx = [2; 41];

    % Format inputs
    [aparcf, out_dir] = format_mask_inputs(aparcf, in_dir, out_dir);
    scan_tag = get_scan_tag(aparcf);

    % Define output files
    outfiles.mask_subcortwm = fullfile(out_dir, append(scan_tag, '_mask-subcortwm.nii'));
    outfiles.ssubcortwm = add_presuf(outfiles.mask_subcortwm, 's');
    outfiles.mask_eroded_subcortwm = fullfile(out_dir, append(scan_tag, '_mask-eroded-subcortwm.nii'));

    % Check if all output files already exist
    if all(structfun(@isfile, outfiles)) && ~overwrite
        log_append(fid, '  * Subcortical white matter mask files exist, will not overwrite');
        return
    end

    % Save the mask
    nii_labels_to_mask(aparcf, mask_idx, outfiles.mask_subcortwm, fid, overwrite);

    % Smooth the mask
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.smooth.data = cellstr(outfiles.mask_subcortwm);
    matlabbatch{1}.spm.spatial.smooth.fwhm = [smooth_by smooth_by smooth_by];
    matlabbatch{1}.spm.spatial.smooth.dtype = spm_type('float32');
    matlabbatch{1}.spm.spatial.smooth.im = 0;  % no implicit masking
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run', matlabbatch);
    log_append(fid, sprintf('  * Saved %s', basename(outfiles.ssubcortwm)));

    % Erode the mask
    outfiles.mask_eroded_subcortwm = fullfile( ...
        out_dir, append(scan_tag, '_mask-eroded-subcortwm.nii') ...
    );
    nii_thresh_to_mask( ...
        outfiles.ssubcortwm, erosion_thresh, Inf, outfiles.mask_eroded_subcortwm, ...
        fid, overwrite ...
    );
end
