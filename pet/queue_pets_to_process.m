function pet_dirs = queue_pets_to_process(log_dir)
    % Return a cell array of processed PET directories
    % for scans that need to be processed
    % ------------------------------------------------------------------
    arguments
        log_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/scans_to_process'
    end

    % Format paths
    log_dir = abspath(log_dir);

    % Print the welcome message
    fprintf('- Loading list of PET scans to process\n');

    % Load the most recent raw_pets table
    filePattern = fullfile(log_dir, 'Raw_PET_Scan_Index_*.csv');
    files = glob_sort_mtime(filePattern);
    fprintf('  * Reading %s\n', files{1});
    raw_pets = readtable(files{1});

    idx = raw_pets.need_to_process == 1;
    pet_dirs = raw_pets.pet_proc_dir(idx);
    fprintf('  * %d/%d PET scans are scheduled for processing\n', sum(idx), length(idx));
end
