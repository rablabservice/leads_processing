function outfiles = process_single_pet(pet_dir, overwrite)
    % Process a single PET scan through all the processing steps.
    %
    % Overview
    % --------
    % 1.  Reset PET origin to the midpoint along each axis.
    % 2.  Coregister and reslice PET to the nu.nii
    % 3.  Calculate reference region means, and save voxelwise PET SUVR
    %     images in native MRI space
    % 4.  Warp PET SUVRs to MNI space using the forward deformation
    %     field from MRI segmentation
    % 5.  Linearly transform PET SUVRs to MNI space using the affine
    %     transform estimated for the nu.nii
    %
    % Usage
    % -----
    % outfiles = process_single_pet(pet_dir, overwrite)
    %
    % Parameters
    % ----------
    % pet_dir : string
    %     The directory containing PET scan data.
    % overwrite : logical
    %     Flag to overwrite existing processed files.
    %
    % Returns
    % -------
    % outfiles : struct
    %     Struct array with paths to each output file
    % ------------------------------------------------------------------
    arguments
        pet_dir {mustBeFolder}
        overwrite logical = false
    end

    % Format paths
    pet_dir = abspath(pet_dir);
    pet_tag = get_scan_tag(pet_dir);
    mri_dir = fullfile(pet_dir, 'mri');
    mri_tag = get_scan_tag(mri_dir);

    % Print the module header
    title = 'PET PROCESSING MODULE';
    subtitle = append('SCAN = ', pet_tag);
    print_title(title, subtitle);

    % Find the raw PET scan
    outfiles.raw_pet = fullfile(pet_dir, append(pet_tag, '.nii'));

    % Initialize SPM jobman and PET parameter defaults
    spm_jobman('initcfg');
    spm('defaults','PET');

    % Reset origin to axis midpoint
    outfiles.raw_pet = reset_origin_axis_midpoint(outfiles.raw_pet, overwrite);

    % Coregister and reslice PET to the nu.nii
    mrifs.nu = fullfile(mri_dir, append(mri_tag, '_nu.nii'));
    outfiles.raw_pet = coreg_reslice_pet_to_mri(outfiles.raw_pet, mrifs.nu, overwrite);

    % Save SUVRs in native MRI space
    outfiles.suvrs = save_pet_suvrs(pet_dir, overwrite, mri_dir);

    % Warp the nu.nii to MNI space using the forward deformation field
    % estimated during segmentation
    mrifs.y = fullfile(mri_dir, append('y_', mri_tag, '_nu.nii'));
    outfiles.wsuvrs = apply_warp_to_mni(outfiles.suvrs, mrifs.y, overwrite);

    % Affine transform the nu.nii to MNI space
    mrifs.atf = fullfile(mri_dir, append('atf_', mri_tag, '_nu.nii'));
    outfiles.apet = apply_affine_to_mni(outfiles.suvrs, mrifs.atf, overwrite);

    % Print the module footer
    print_footer('PET processing module complete');
end
