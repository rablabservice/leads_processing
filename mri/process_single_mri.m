function process_single_mri(mri_dir, skip_freesurfer, segment_brainstem, overwrite, verbose)
    % Run a single MRI scan through all the processing steps.
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
        skip_freesurfer logical = false
        segment_brainstem logical = true
        overwrite logical = false
        verbose logical = true
    end

    % Format paths
    mri_dir = abspath(mri_dir);

    % Print the scan header
    if verbose
        scan_tag = get_scan_tag(mri_dir);
        fprintf('\n%s\n', scan_tag);
        fprintf('%s\n', repmat('-', 1, length(scan_tag)));
    end

    % Run FreeSurfer
    if ~skip_freesurfer
        raw_mrif = fullfile(mri_dir, 'raw', 'mri.nii.gz');
        process_mri_freesurfer(raw_mrif, mri_dir, segment_brainstem, overwrite, verbose);
    end

    % Run post-FreeSurfer processing
    process_mri_post_freesurfer(mri_dir, segment_brainstem, overwrite, verbose);
end
