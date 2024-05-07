function outfiles = process_single_mri(mri_dir, segment_brainstem, overwrite)
    % Run a single MRI scan through all the processing steps.
    %
    % Usage
    % -----
    % outfiles = process_single_mri(mri_dir)
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
        segment_brainstem logical = true
        overwrite logical = false
    end

    % Format paths
    mri_dir = abspath(mri_dir);
    scan_tag = get_scan_tag(mri_dir);

    % Print the module header
    title = 'MRI PROCESSING MODULE';
    subtitle = append('SCAN = ', scan_tag);
    print_title(title, subtitle);

    % % Run FreeSurfer
    % raw_mrif = fullfile(mri_dir, 'raw', 'mri.nii.gz');
    % process_mri_freesurfer(raw_mrif, mri_dir, segment_brainstem, overwrite);

    % Run post-FreeSurfer processing
    outfiles = process_mri_post_freesurfer(mri_dir, overwrite);
end
