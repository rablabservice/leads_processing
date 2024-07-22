function outfiles = process_mri_post_freesurfer(mri_dir, fid, overwrite)
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
    % 5.  Segment nu.nii and save deformation fields between native MRI
    %     and MNI space
    % 6.  Save reference region and target ROI masks in native MRI space
    % 7.  Warp nu.nii to MNI space
    % 8.  Use nu.nii to calculate and save the full affine transform
    %     from native MRI to MNI space
    % 9.  Apply affine transformation to the native space nu.nii
    %
    % Parameters
    % ----------
    % mri_dir : char or str
    %     Path to the subject's processed MRI directory
    % fid : int, optional
    %     File identifier for logging (default is 1 for stdout)
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
        fid {mustBeNumeric} = 1
        overwrite logical = false
    end

    % Format inputs
    mri_dir = abspath(mri_dir);

    % If mri_dir is not the baseline MRI and the baseline MRI has not
    % already been processed, don't process this MRI
    subj_dir = fileparts(mri_dir);
    baseline_mri_dir = fileparts(get_baseline_mri(subj_dir));
    if ~strcmp(mri_dir, baseline_mri_dir) && ~processed_mri_files_exist(baseline_mri_dir)
        log_append(fid, append('- Cannot proceed with post-FreeSurfer MRI processing, because the baseline MRI has not been fully processed'));
        outfiles = struct([]);
        return
    end

    % Initialize SPM jobman and PET parameter defaults
    spm_jobman('initcfg');
    spm('defaults','PET');

    % Check if FreeSurfer files have already been converted to nifti
    outfiles = get_freesurfer_files(mri_dir);
    if all(structfun(@isfile, outfiles)) && ~overwrite
        log_append(fid, sprintf([ ...
            '- FreeSurfer .mgz files have already been converted to .nii; will skip\n' ...
            '            conversion, recentering, and coregistering these files to the baseline\n' ...
            '            MRI as these steps have presumably already been completed. Rerun with\n' ...
            '            overwrite=true to redo these steps' ...
        ]));
        skip_mri_convert_recenter_coreg = true;
    else
        skip_mri_convert_recenter_coreg = false;
    end


    % Copy FreeSurfer files to the subject's processed MRI directory
    % and convert them from .mgz to .nii
    if ~skip_mri_convert_recenter_coreg
        outfiles = copy_convert_freesurfer(mri_dir, fid, overwrite);
    end

    % Reset origin of the nu.nii to its center-of-mass, then coregister
    % to the OldNorm T1 template. Apply transforms to aparc+aseg and
    % brainstem seg files. This step overwrites affine transforms in the
    % input image headers but does not affect their data arrays
    if ~skip_mri_convert_recenter_coreg
        reset_origin_mri_com(outfiles, fid);
    end

    % Coregister nu.nii to subject's baseline nu.nii. Apply transforms
    % to aparc+aseg and brainstem seg files. Overwrites affines in the
    % input image headers but does not affect their data arrays
    if ~skip_mri_convert_recenter_coreg
        coreg_mri_to_baseline(outfiles, fid);
    end

    % Run SUIT to get the cerebellar atlas in native MRI space
    outfiles = catstruct(outfiles, run_suit(mri_dir, fid, overwrite));

    % Segment the nu.nii and save forward and inverse deformation fields
    % for later use in warping images between native MRI and MNI space
    outfiles = catstruct(outfiles, segment_mri(outfiles.nu, fid, overwrite));

    % Save reference region and target ROI mask files
    outfiles = catstruct(outfiles, save_roi_masks(mri_dir, fid, overwrite));

    % Warp the nu.nii to MNI space using the forward deformation field
    % estimated during segmentation
    outfiles.wnu = apply_warp_to_mni( ...
        outfiles.nu, outfiles.y, fid, overwrite ...
    );

    % Calculate affine transform from native MRI to MNI space
    outfiles.atf = estimate_mri_affine_to_mni(outfiles.nu, fid, overwrite);

    % Affine transform the nu.nii to MNI space
    outfiles.anu = apply_affine_to_mni( ...
        outfiles.nu, outfiles.atf, fid, overwrite ...
    );
end
