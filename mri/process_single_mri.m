function process_single_mri(scan_tag, data_dir, overwrite, verbose)
    % Run a single MRI scan through all the processing steps.
    % ------------------------------------------------------------------

    % Get the scan tag
    [subj, scan_type, scan_date] = parse_scan_tag(scan_tag);


    % Find the scan in raw


    % Look for the scan in processed


    % Run FreeSurfer processing


    % Run post-FreeSurfer processing

end
