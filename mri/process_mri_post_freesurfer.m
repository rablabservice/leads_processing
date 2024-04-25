function process_mri_post_freesurfer(mri_dir, segment_brainstem, overwrite, verbose)
    % Complete post-FreeSurfer MRI processing, mainly SPM12 operations
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeText}
        segment_brainstem logical = true
        overwrite logical = false
        verbose logical = true
    end

    % Format inputs
    mri_dir = abspath(mri_dir);
    subj_dir = fileparts(mri_dir);

    % Copy FreeSurfer files to the subject's processed mri directory
    % and convert them from .mgz to .nii
    mri_files = cell(1, 2 + segment_brainstem);
    [mri_files{:}] = copy_convert_freesurfer(mri_dir, segment_brainstem, overwrite, verbose);

    % Reset origin to center-of-mass (note this step overwrites the
    % input image files)
    % mri_files = mri_reset_origin_com(mri_files, verbose);

    % If MRI is not baseline, coregister it to the baseline MRI
    baseline_nuf = get_baseline_mri(subj_dir);
    if ~strcmp(mri_files{1}, baseline_nuf)
        error('MRI is not baseline, coregistration to baseline MRI not yet implemented')
        baseline_mri_dir = get_baseline_mri_dir(mri_dir);
        mri_files = coreg_mri_to_baseline(mri_files, baseline_mri_dir);
    end

    % % Segment MRI, save forward and inverse deformation fields, and warp
    % % the nu.nii to MNI space
    % segment_mri_and_warp_to_mni(mri_dir, overwrite, verbose);

    % % Calculate affine transform from
    % affine_mri_to_mni(mri_dir, overwrite, verbose);

    % Save out mask files used as reference regions or target ROIs.
    save_aparc_roi_masks(mri_dir, overwrite, verbose);
end



