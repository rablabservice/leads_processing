function outfiles = process_mri_post_freesurfer(mri_dir, overwrite)
    % Complete post-FreeSurfer MRI processing, mainly SPM12 operations
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeText}
        overwrite logical = false
    end

    % Format inputs
    mri_dir = abspath(mri_dir);
    subj_dir = fileparts(mri_dir);

    % Initialize SPM jobman and PET parameter defaults
    spm_jobman('initcfg');
    spm('defaults','PET');

    % Copy FreeSurfer files to the subject's processed mri directory
    % and convert them from .mgz to .nii
    outfiles = copy_convert_freesurfer(mri_dir, overwrite);

    % Reset origin of the nu.nii to its center-of-mass, then coregister
    % to the OldNorm T1 template. Apply transforms to aparc+aseg and
    % brainstem seg files. This step overwrites affine transforms in the
    % input image headers but does not affect their data arrays
    reset_origin_mri_com(outfiles, '', false);

    % Coregister nu.nii to subject's baseline nu.nii. Apply transforms
    % to aparc+aseg and brainstem seg files. Overwrites affines in the
    % input image headers but does not affect their data arrays
    coreg_mri_to_baseline(outfiles);

    % % Segment the nu.nii and save forward and reverse deformation fields
    % % for later use in warping images from native MRI to MNI space
    outfiles = catstruct(outfiles, segment_mri(outfiles.nu, overwrite));

    % Warp nu to MNI space
    interp = 4;
    vox = [1.5 1.5 1.5];
    prefix = 'w';
    bb = [Inf Inf Inf; Inf Inf Inf];
    outfiles.wnu = apply_warp_to_mni(outfiles.nu, outfiles.y, interp, vox, prefix, bb, overwrite);

    % % Calculate affine transform from MRI native space to MNI
    % estimate_mri_affine_to_mni(mri_dir, overwrite);

    % Run SUIT to get the cerebellar atlas in native MRI space
    outfiles = catstruct(outfiles, run_suit(mri_dir, overwrite));

    % Save reference region and target ROI mask files
    outfiles = catstruct(outfiles, save_roi_masks(mri_dir, overwrite));
end
