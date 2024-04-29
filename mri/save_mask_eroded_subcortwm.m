function outfiles = save_mask_eroded_subcortwm( ...
    aparcf, smooth_by, erosion_thresh, in_dir, out_dir, overwrite, verbose ...
)
    % Load the aparc+aseg and save eroded subcortical white matter mask
    %
    % Usage
    % -----
    % >> save_mask_eroded_subcortwm(aparcf)
    %
    % Parameters
    % ----------
    % aparcf : char or str array
    %   Path to the aparc+aseg.nii file
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
    % overwrite : logical, optional
    %   If true, overwrite existing file
    % verbose : logical, optional
    %   If true, print diagnostic information
    %
    % Output
    % ------

    %
    % Files created
    % -------------
    % - <mri_dir>/<scan_tag>_mask-eroded-subcortwm.nii
    % - <mri_dir>/<scan_tag>_mask-subcortwm.nii
    % - <mri_dir>/s<scan_tag>_mask-subcortwm.nii
    % ------------------------------------------------------------------
    arguments
        aparcf {mustBeText} = ''
        smooth_by {mustBeNumeric} = 8
        erosion_thresh {mustBeNumeric} = 0.7
        in_dir {mustBeText} = ''
        out_dir {mustBeText} = ''
        overwrite logical = false
        verbose logical = true
    end

    % Define aparc indices for the subcortical white matter
    mask_idx = [2; 41];

    % Format inputs
    [aparcf, out_dir] = format_mask_inputs(aparcf, in_dir, out_dir);
    scan_tag = get_scan_tag(aparcf);

    % Save the mask
    outfiles.subcortwm = fullfile(out_dir, append(scan_tag, '_mask-subcortwm.nii'));
    nii_labels_to_mask(aparcf, mask_idx, outfiles.subcortwm, overwrite, verbose);

    % Smooth the mask
    outfiles.ssubcortwm = fullfile( ...
        out_dir, append('s', scan_tag, '_mask-subcortwm.nii') ...
    );
    clear matlabbatch;
    matlabbatch{2}.spm.spatial.smooth.data = cellstr(outfiles.subcortwm);
    matlabbatch{2}.spm.spatial.smooth.fwhm = [smooth_by smooth_by smooth_by];
    matlabbatch{2}.spm.spatial.smooth.dtype = 0;  % same as input
    matlabbatch{2}.spm.spatial.smooth.im = 0;  % no implicit masking
    matlabbatch{2}.spm.spatial.smooth.prefix = 's';
    spm('defaults', 'PET');
    spm_jobman('run', matlabbatch);

    % Erode the mask
    outfiles.eroded_subcortwm = fullfile( ...
        out_dir, append(scan_tag, '_mask-eroded-subcortwm.nii') ...
    );
    nii_thresh_to_mask( ...
        outfiles.ssubcortwm, erosion_thresh, Inf, outfiles.eroded_subcortwm, ...
        overwrite, verbose ...
    );
end
