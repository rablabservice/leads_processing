function mri_dirs = select_mris_to_process(log_dir, verbose)
    % Return an array of MRIs to process by comparing raw to processed.
    % ------------------------------------------------------------------
    arguments
        log_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/log'
        verbose logical = true
    end

    % Format paths
    log_dir = abspath(log_dir);

    % Load the most recent raw_scans table
    filePattern = fullfile(log_dir, 'raw_pet_scans_*.csv');
    files = glob_sort_mtime(filePattern);
    raw_scans = readtable(files{1});

    % Get all unique MRI directories
    mri_dirs = unique(raw_scans.mri_dir);

    if verbose
        fprintf('- Found %d unique MRI directories\n', length(mri_dirs));
    end
end
