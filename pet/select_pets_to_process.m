function pet_dirs = select_pets_to_process(log_dir)
    % Return an array of PET scans to process
    % ------------------------------------------------------------------
    arguments
        log_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/scans_to_process'
    end

    % Format paths
    log_dir = abspath(log_dir);

    % Load the most recent raw_pets table
    filePattern = fullfile(log_dir, 'Raw_PET_Scan_Index_*.csv');
    files = glob_sort_mtime(filePattern);
    raw_pets = readtable(files{1});

    idx = raw_pets.need_to_process == 1;
    pet_dirs = raw_pets.pet_proc_dir(idx);

    fprintf('- Selected %d PET scans to process\n', sum(idx));
end
