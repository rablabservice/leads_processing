function outfiles = save_roi_masks(mri_dir, overwrite, verbose)
    % Save out mask files used as reference regions or target ROIs.
    %
    % Parameters
    % ----------
    % mri_dir : char or str array
    %   The directory that contains the processed MRI data
    % overwrite : logical, optional
    %   If true, overwrite existing files
    % verbose : logical, optional
    %   If true, print diagnostic information
    %
    % Files created
    % -------------
    % - <mri_dir>/<scan_tag>_cbl-suit.nii
    % - <mri_dir>/<scan_tag>_mask-amyloid-cortical-summary.nii
    % - <mri_dir>/<scan_tag>_mask-brainstem.nii
    % - <mri_dir>/<scan_tag>_mask-eroded-subcortwm.nii
    % - <mri_dir>/<scan_tag>_mask-infcblgm.nii
    % - <mri_dir>/<scan_tag>_mask-pons.nii
    % - <mri_dir>/<scan_tag>_mask-wcbl.nii
    % - <mri_dir>/<scan_tag>_mask-subcortwm.nii
    % - <mri_dir>/s<scan_tag>_mask-subcortwm.nii
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
        overwrite logical = false
        verbose logical = true
    end

    % Format paths
    mri_dir = abspath(mri_dir);
    [~, ~, aparcf, brainstemf] = get_freesurfer_files(mri_dir);

    % Check that atlas files exist
    save_aparc_masks = exist(aparcf, 'file');
    save_brainstem_masks = exist(brainstemf, 'file')
    if ~save_aparc_masks && ~save_brainstem_masks
        if verbose
            fprintf('- No aparc or brainstem sublabel files found in %s\n', mri_dir);
        end
        return
    end

    % Initialize output
    outfiles = struct([]);

    % Call the individual ROI creation functions
    if verbose
        fprintf('- Saving ROI masks\n');
    end
    if save_aparc_masks
        outfiles.wcbl = save_mask_wcbl(aparcf, '', '', overwrite, verbose);
        outfiles.brainstem = save_mask_brainstem(aparcf, '', '', overwrite, verbose);
        outfiles.cortical_summary = save_mask_amyloid_cortical_summary(aparcf, '', '', overwrite, verbose);
        outfiles = catstruct(outfiles, save_mask_eroded_subcortwm(aparcf, 8, 0.7, '', '', overwrite, verbose));
        outfiles = catstruct(outfiles, save_mask_infcblgm(aparcf, '', '', overwrite, verbose));
    end
    if save_brainstem_masks
        outfiles.pons = save_mask_pons(brainstemf, '', '', overwrite, verbose);
    end
end
