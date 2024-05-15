function process_pets(overwrite, log_dir)
    % High-level function to select and process PET scans
    %
    % Parameters
    % ----------
    % overwrite : logical
    %     If true, overwrite existing processed data
    % data_dir : char or str
    %     The directory that contains raw (unprocessed) PET data in
    %     <data_dir>/raw and processed data in <data_dir>/processed
    % log_dir : char or str
    %     The directory that stores log files
    % ------------------------------------------------------------------
    arguments
        overwrite logical = false
        log_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/scans_to_process'
    end

    % Format paths
    log_dir = abspath(log_dir);

    % Select MRIs to process
    pet_dirs = select_pets_to_process(log_dir);

    % Process MRIs in parallel
    parfor ii = 1:length(pet_dirs)
        process_single_pet(pet_dirs(ii), overwrite);
    end
end
