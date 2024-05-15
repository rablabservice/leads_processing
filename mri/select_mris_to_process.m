function [raw_mrifs, mri_dirs] = select_mris_to_process(log_dir)
    % Return an array of MRIs to process
    % ------------------------------------------------------------------
    arguments
        log_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/scans_to_process'
    end

    % Format paths
    log_dir = abspath(log_dir);

    % Load the most recent raw_MRIs table
    filePattern = fullfile(log_dir, 'Raw_MRI_Scan_Index_*.csv');
    files = glob_sort_mtime(filePattern);
    raw_mris = readtable(files{1});

    idx = raw_mris.need_to_process == 1;
    raw_mrifs = raw_mris.mri_raw_niif(idx);
    mri_dirs = raw_mris.mri_proc_dir(idx);

    fprintf('- Selected %d MRIs to process\n', sum(idx));
end
