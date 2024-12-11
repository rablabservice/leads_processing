function process_pets(pet_dirs, scans_to_process_dir, overwrite, run_qc)
    % Process all PET scans that are scheduled for processing
    % in the latest raw_PET_index file
    %
    % Parameters
    % ----------
    % pet_dirs : cell
    %     List of PET directories to process. This variable overrides
    %     the default behavior of selecting scans to process based on
    %     the latest raw_PET_index file in scans_to_process_dir
    % scans_to_process_dir : char or str
    %     The directory that stores raw_PET_index files
    % overwrite : logical
    %     If true, overwrite existing processed data
    % run_qc : logical, optional
    %     If true, create QC image and add new QC eval file. Default is
    %     true
    % ------------------------------------------------------------------
    arguments
        pet_dirs = {}
        scans_to_process_dir {mustBeText} = '/mnt/coredata/processing/leads/metadata/scans_to_process'
        overwrite logical = false
        run_qc logical = true
    end

    % Format paths
    scans_to_process_dir = abspath(scans_to_process_dir);

    % Load the list of PET scans to process
    fprintf('Checking list of submitted scans to process...\n');
    if isempty(pet_dirs)
        pet_dirs = queue_pets_to_process(scans_to_process_dir);
    else
        pet_dirs = abspath(cellvec(pet_dirs));
        pet_dirs = pet_dirs(cellfun(@isfolder, pet_dirs));
    end

    % If overwrite is false, remove scans from the list that have
    % already been processed
    if ~overwrite
        pet_dirs = pet_dirs(cellfun(@(x) ~processed_pet_files_exist(x), pet_dirs));
    end

    % Are there are any scans to process?
    if isempty(pet_dirs)
        fprintf('No PET scans to process\n');
        return;
    elseif length(pet_dirs) == 1
        fprintf('Preparing to process 1 PET scan\n');
        cellfun(@(x) fprintf('  %s\n', x), pet_dirs);
    else
        fprintf('Preparing to process %d PET scans\n', length(pet_dirs));
        cellfun(@(x) fprintf('  %s\n', x), pet_dirs);
    end

    % Process one scan
    if length(pet_dirs) == 1
        process_single_pet(pet_dirs{1}, overwrite, run_qc);
    % Process multiple scans in parallel
    else
        % Start a parallel pool
        n_workers = min(length(pet_dirs), maxNumCompThreads);
        poolobj = parpool(n_workers);

        % Assign a worker to each scan
        parfor ii = 1:length(pet_dirs)
            try
                process_single_pet(pet_dirs{ii}, overwrite, run_qc);
            catch ME
                warning( ...
                    '\n\n\nERROR processing %s: %s\n\n\n', pet_dirs{ii}, ...
                    getReport(ME, 'extended', 'hyperlinks', 'off') ...
                );
            end
        end

        % Close the parallel pool
        if ~isempty(poolobj)
            delete(poolobj);
        end
    end
end
