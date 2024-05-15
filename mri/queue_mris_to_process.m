function [raw_mrifs, mri_dirs] = queue_mris_to_process(log_dir)
    % Return cell arrays of raw MRI files and processed MRI directories
    % for scans that need to be processed
    % ------------------------------------------------------------------
    arguments
        log_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/scans_to_process'
    end

    % Format paths
    log_dir = abspath(log_dir);

    % Print the welcome message
    fprintf('- Loading list of MRIs to process\n');

    % Load the most recent raw_MRIs table
    filePattern = fullfile(log_dir, 'Raw_MRI_Scan_Index_*.csv');
    files = glob_sort_mtime(filePattern);
    fprintf('  * Reading %s\n', files{1});
    raw_mris = readtable(files{1});

    idx = raw_mris.need_to_process == 1;
    raw_mrifs = raw_mris.mri_raw_niif(idx);
    mri_dirs = raw_mris.mri_proc_dir(idx);
    fprintf('  * %d/%d MRIs are scheduled for processing\n', sum(idx), length(idx));
end
