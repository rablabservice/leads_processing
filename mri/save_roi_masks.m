function outfiles = save_roi_masks(mri_dir, aparcf, bstemf, overwrite)
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
    % aparcf : char or str array, optional
    %   Path to the aparc+aseg .nii file
    % bstemf : char or str array, optional
    %   Path to the brainstem sublabels .nii file
    % overwrite : logical, optional
    %   If true, overwrite existing files
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
        aparcf = ''
        bstemf = ''
        overwrite logical = false
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
    iyf = fullfile(mri_dir, append('iy_', scan_tag, '_nu.nii'));

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
        outfiles = catstruct(outfiles, save_mask_wcbl(aparcf, '', '', overwrite));
        outfiles = catstruct(outfiles, save_mask_brainstem(aparcf, '', '', overwrite));
        outfiles = catstruct(outfiles, save_mask_amyloid_cortical_summary(aparcf, '', '', overwrite));
        outfiles = catstruct(outfiles, save_mask_eroded_subcortwm(aparcf, 8, 0.7, '', '', overwrite));
        if isfile(iyf)
            outfiles = catstruct(outfiles, save_mask_infcblgm(aparcf, iyf, '', '', overwrite));
        end
    end
    if save_brainstem_masks
        outfiles.pons = save_mask_pons(bstemf, '', '', overwrite);
    end
end
