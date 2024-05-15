function process_pets(overwrite, log_dir)
    % Process all PET scans that are scheduled for processing
    % in the latest log file
    %
    % Parameters
    % ----------
    % overwrite : logical
    %     If true, overwrite existing processed data
    % log_dir : char or str
    %     The directory that stores log files
    % ------------------------------------------------------------------
    arguments
        overwrite logical = false
        log_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/scans_to_process'
    end

    % Format paths
    log_dir = abspath(log_dir);

    % Load the list of PET scans to process
    pet_dirs = queue_pets_to_process(log_dir);

    % Process PET scans in parallel
    parfor ii = 1:length(pet_dirs)
        process_single_pet(pet_dirs(ii), overwrite);
    end
end
