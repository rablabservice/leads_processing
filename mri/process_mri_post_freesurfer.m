function process_mri_post_freesurfer(mri_dir, segment_brainstem, overwrite, verbose)
    % Complete post-FreeSurfer MRI processing, mainly SPM12 operations
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeText}
        segment_brainstem logical = true
        overwrite logical = false
        verbose logical = true
    end

    % Format inputs; initialize SPM jobman and PET parameter defaults
    mri_dir = abspath(mri_dir);
    subj_dir = fileparts(mri_dir);
    spm_jobman('initcfg');
    spm('defaults','PET');

    % Copy FreeSurfer files to the subject's processed mri directory
    % and convert them from .mgz to .nii
    mri_files = cell(1, 2 + segment_brainstem);
    [mri_files{:}] = copy_convert_freesurfer( ...
        mri_dir, segment_brainstem, overwrite, verbose ...
    );

    % Reset origin of the nu.nii to its center-of-mass, then coregister
    % it to the OldNorm T1 template. Apply transforms to aparc+aseg and
    % brainstem seg files. This step overwrites affine transforms in the
    % input image headers but does not affect their data arrays
    mri_files = reset_origin_mri_com(mri_files, verbose);

    % % Coregister nu.nii to subject's baseline nu.nii. Apply transforms
    % % to aparc+aseg and brainstem seg files. Overwrites affines in the
    % % input image headers but does not affect their data arrays
    % mri_files = coreg_mri_to_baseline(mri_files);

    % % Segment the nu.nii and save forward and reverse deformation fields
    % % for later use in warping images from subject MRI to MNI space
    % segment_mri(mri_dir, overwrite, verbose);

    % % Warp MRI files to MNI space
    % warp_to_mni(nuf, mri_dir, overwrite, verbose);

    % % Calculate affine transform from MRI native space to MNI
    % affine_mri_to_mni(mri_dir, overwrite, verbose);

    % Save reference region and target ROI mask files
    save_roi_masks(mri_dir, overwrite, verbose);
end
