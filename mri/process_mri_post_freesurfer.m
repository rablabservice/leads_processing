function outfiles = process_mri_post_freesurfer(mri_dir, overwrite)
    % Complete post-FreeSurfer MRI processing
    %
    % Overview
    % --------
    % 1.  Copy and convert FreeSurfer files to NIfTI
    % 2.  Reset FreeSurfer file origins to the nu.nii center-of-mass,
    %     then coregister nu.nii to MNI space
    % 3.  Coregister nu.nii to the baseline nu.nii (earliest MRI for a
    %     subject). This step is skipped if nu.nii is at baseline.
    %     FreeSurfer files are now in native MRI space
    % 4.  Run SUIT and inverse warp cerebellar subregions to native MRI
    %     space
    % 5.  Save reference region and target ROI masks in native MRI space
    % 6.  Segment nu.nii and save deformation fields between native MRI
    %     and MNI space
    % 7.  Warp nu.nii to MNI space
    % 8.  Use nu.nii to calculate and save the full affine transform
    %     from native MRI to MNI space
    % 9.  Apply affine transformation to the native space nu.nii
    %
    % Parameters
    % ----------
    % mri_dir : char or str
    %     Path to the subject's processed MRI directory
    % overwrite : logical, optional
    %     If true, overwrite existing files (default is false)
    %
    % Returns
    % -------
    % outfiles : struct
    %     Struct array with paths to each output file
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeText}
        overwrite logical = false
    end

    % Format inputs
    mri_dir = abspath(mri_dir);

    % Initialize SPM jobman and PET parameter defaults
    spm_jobman('initcfg');
    spm('defaults','PET');

    % Copy FreeSurfer files to the subject's processed MRI directory
    % and convert them from .mgz to .nii
    outfiles = copy_convert_freesurfer(mri_dir, overwrite);

    % Reset origin of the nu.nii to its center-of-mass, then coregister
    % to the OldNorm T1 template. Apply transforms to aparc+aseg and
    % brainstem seg files. This step overwrites affine transforms in the
    % input image headers but does not affect their data arrays
    reset_origin_mri_com(outfiles);

    % Coregister nu.nii to subject's baseline nu.nii. Apply transforms
    % to aparc+aseg and brainstem seg files. Overwrites affines in the
    % input image headers but does not affect their data arrays
    coreg_mri_to_baseline(outfiles);

    % Run SUIT to get the cerebellar atlas in native MRI space
    outfiles = catstruct(outfiles, run_suit(mri_dir, overwrite));

    % Save reference region and target ROI mask files
    outfiles = catstruct(outfiles, save_roi_masks(mri_dir, overwrite));

    % Segment the nu.nii and save forward and inverse deformation fields
    % for later use in warping images between native MRI and MNI space
    outfiles = catstruct(outfiles, segment_mri(outfiles.nu, overwrite));

    % Warp the nu.nii to MNI space using the forward deformation field
    % estimated during segmentation
    outfiles.wnu = apply_warp_to_mni( ...
        outfiles.nu, outfiles.y, overwrite ...
    );

    % Calculate affine transform from native MRI to MNI space
    outfiles.atf = estimate_mri_affine_to_mni(outfiles.nu, overwrite);

    % Affine transform the nu.nii to MNI space
    outfiles.anu = apply_affine_to_mni( ...
        outfiles.nu, outfiles.atf, overwrite ...
    );
end
