function outfiles = process_single_mri(raw_mrif, mri_dir, overwrite, segment_brainstem)
    % Run a single MRI scan through all the processing steps.
    %
    % Overview
    % --------
    % 1.  Run recon-all on the raw MRI file
    % 2.  Optionally segment the brainstem using segmentBS.sh
    % 3.  Copy and convert FreeSurfer files to NIfTI
    % 4.  Reset FreeSurfer file origins to the nu.nii center-of-mass,
    %     then coregister nu.nii to MNI space
    % 5.  Coregister nu.nii to the baseline nu.nii (earliest MRI for a
    %     subject). This step is skipped if nu.nii is at baseline.
    %     FreeSurfer files are now in native MRI space
    % 6.  Run SUIT and inverse warp cerebellar subregions to native MRI
    %     space
    % 7.  Save reference region and target ROI masks in native MRI space
    % 8.  Segment nu.nii and save deformation fields between native MRI
    %     and MNI space
    % 9.  Warp nu.nii to MNI space
    % 10. Use nu.nii to calculate and save the full affine transform
    %     from native MRI to MNI space
    % 11. Apply affine transformation to the native space nu.nii
    %
    % Usage
    % -----
    % outfiles = process_single_mri(mri_dir)
    %
    % Parameters
    % ----------
    % raw_mrif : char or str array
    %     Full path to the raw MRI nifti
    % mri_dir : char or str array
    %     Path to the processed MRI directory
    % overwrite : logical, optional
    %     If true, overwrite existing files. Default is false
    % segment_brainstem : logical, optional
    %     If true, segment the brainstem using segmentBS.sh. Default is
    %     true
    % ------------------------------------------------------------------
    arguments
        raw_mrif {mustBeFile}
        mri_dir {mustBeFolder}
        overwrite logical = false
        segment_brainstem logical = true
    end

    % Format paths
    mri_dir = abspath(mri_dir);
    scan_tag = get_scan_tag(mri_dir);

    % Print the module header
    title = 'MRI PROCESSING MODULE';
    subtitle = append('SCAN : ', scan_tag);
    print_title(title, subtitle);

    % Run FreeSurfer
    process_mri_freesurfer(raw_mrif, mri_dir, overwrite, segment_brainstem);

    % Run post-FreeSurfer processing
    outfiles = process_mri_post_freesurfer(mri_dir, overwrite);

    % Print the module footer
    print_footer('MRI processing module complete');
end
