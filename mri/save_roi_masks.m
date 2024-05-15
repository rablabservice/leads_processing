function outfiles = save_roi_masks(mri_dir, overwrite, aparcf, bstemf)
    % Save out mask files used as reference regions or target ROIs.
    %
    % If mri_dir is passed, look up the aparc and brainstem files in
    % mri_dir following expected nomenclature. If mri_dir is empty and
    % aparcf and/or bstemf are passed, use these paths directly.
    %
    % Parameters
    % ----------
    % mri_dir : char or str array
    %   The directory that contains the processed MRI data
    % overwrite : logical, optional
    %   If true, overwrite existing files
    % aparcf : char or str array, optional
    %   Path to the aparc+aseg .nii file
    % bstemf : char or str array, optional
    %   Path to the brainstem sublabels .nii file
    %
    % Files created
    % -------------
    % - <mri_dir>/<scan_tag>_cbl-suit.nii
    % - <mri_dir>/<scan_tag>_mask-amyloid-cortical-summary.nii
    % - <mri_dir>/<scan_tag>_mask-brainstem.nii
    % - <mri_dir>/<scan_tag>_mask-cblgm.nii
    % - <mri_dir>/<scan_tag>_mask-eroded-subcortwm.nii
    % - <mri_dir>/<scan_tag>_mask-infcblgm.nii
    % - <mri_dir>/<scan_tag>_mask-pons.nii
    % - <mri_dir>/<scan_tag>_mask-wcbl.nii
    % - <mri_dir>/<scan_tag>_mask-subcortwm.nii
    % - <mri_dir>/s<scan_tag>_mask-subcortwm.nii
    % ------------------------------------------------------------------
    arguments
        mri_dir = ''
        overwrite logical = false
        aparcf = ''
        bstemf = ''
    end

    % Format parameters
    if ~isempty(mri_dir)
        mrifs = get_freesurfer_files(mri_dir);
        aparcf = mrifs.aparc;
        bstemf = mrifs.bstem;
    else
        aparcf = abspath(aparcf);
        bstemf = abspath(bstemf);
    end
    scan_tag = get_scan_tag(mri_dir);
    suitf = fullfile(mri_dir, append(scan_tag, '_cbl-suit.nii'));

    % Check that atlas files exist
    save_aparc_masks = isfile(aparcf);
    save_brainstem_masks = isfile(bstemf);
    if ~save_aparc_masks && ~save_brainstem_masks
        fprintf('- No aparc or brainstem sublabel files found in %s\n', mri_dir);
        return
    end

    % Initialize output
    outfiles = struct();

    % Call the individual ROI creation functions
    fprintf('- Saving ROI masks\n');
    if save_aparc_masks
        outfiles = catstruct(outfiles, save_mask_wcbl(aparcf, overwrite));
        outfiles = catstruct(outfiles, save_mask_brainstem(aparcf, overwrite));
        outfiles = catstruct(outfiles, save_mask_amyloid_cortical_summary(aparcf, overwrite));
        outfiles = catstruct(outfiles, save_mask_eroded_subcortwm(aparcf, overwrite));
        if isfile(suitf)
            outfiles = catstruct(outfiles, save_mask_infcblgm(aparcf, suitf, overwrite));
        end
    end
    if save_brainstem_masks
        outfiles.pons = save_mask_pons(bstemf, overwrite);
    end
end
