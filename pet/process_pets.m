function process_pets(overwrite, scans_to_process_dir)
    % Process all PET scans that are scheduled for processing
    % in the latest log file
    %
    % Parameters
    % ----------
    % overwrite : logical
    %     If true, overwrite existing processed data
    % scans_to_process_dir : char or str
    %     The directory that stores log files
    % ------------------------------------------------------------------
    arguments
        overwrite logical = false
        scans_to_process_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/scans_to_process'
    end

    % Format paths
    scans_to_process_dir = abspath(scans_to_process_dir);

    % Load the list of PET scans to process
    pet_dirs = queue_pets_to_process(scans_to_process_dir);

    % Process PET scans in parallel
    parfor ii = 1:length(pet_dirs)
        try
            process_single_pet(pet_dirs{ii}, overwrite);
        catch ME
            warning('\n\n\nERROR processing %s: %s\n\n\n', pet_dirs{ii}, ME.message);
        end
    end
end
