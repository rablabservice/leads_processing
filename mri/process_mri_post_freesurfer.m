function process_mri_post_freesurfer(mri_dir, segment_brainstem, overwrite, verbose)
    % High-level function to process MRI after FreeSurfer has been run
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeText}
        segment_brainstem logical = true
        overwrite logical = false
        verbose logical = true
    end

    % Format inputs
    mri_dir = abspath(mri_dir);
    subj_dir = dirname(mri_dir);
    scan_tag = get_scan_tag(mri_dir);

    % Copy FreeSurfer files to the subject's processed mri directory
    % and convert them from .mgz to .nii
    copy_convert_freesurfer(mri_dir, segment_brainstem, overwrite, verbose);
    if segment_brainstem
        [~, nuf, aparcf, brainstemf] = get_freesurfer_files(mri_dir, 'nii', segment_brainstem);
        mri_files = {nuf, aparcf, brainstemf};
    else
        [~, nuf, aparcf] = get_freesurfer_files(mri_dir, 'nii', segment_brainstem);
        mri_files = {nuf, aparcf};
    end

    % Reset origin to center-of-mass
    mri_reset_origin_com(mri_files);  % overwrites the input images

    % If MRI is not baseline, coregister it to the baseline MRI
    baseline_nuf = get_baseline_mri(subj_dir);

    if ~is_baseline_mri(mri_dir)
        baseline_mri_dir = get_baseline_mri_dir(mri_dir);
        mri_coregister(mri_files, baseline_mri_dir);
    end

    % Save out mask files used as reference regions or target ROIs.


    % Segment MRI, save forward and inverse deformation fields, and warp
    % the nu.nii to MNI space


    % Calculate affine transform from



end



