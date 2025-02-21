function process_wmaps(...
    run_dirs, ...
    scans_to_process_dir, ...
    overwrite ...
)
 % Schedule scans to process W-maps
    %
    % Parameters
    % ----------
    % run_dirs : cell
    %     List of MRI ad/or PET directories to process. This variable overrides
    %     the default behavior of selecting scans to process based on
    %     the latest files in scans_to_process_dir
    % scans_to_process_dir : char or str
    %     The directory that stores scheduled MRI and PET scans list
    % overwrite : logical
    %     If true, overwrite existing W-maps
    % ------------------------------------------------------------------
    arguments
        run_dirs = {}
        scans_to_process_dir {mustBeText} = '/mnt/coredata/processing/leads/metadata/scans_to_process'
        overwrite logical = false
    end

    % Format paths
    scans_to_process_dir = abspath(scans_to_process_dir);

    % Load the list of MRIs to process
    fprintf('Checking list of submitted scans to process...\n');
    if isempty(run_dirs)
        mri_dirs = {};
        pet_dirs = {};

        % Schedule MRIs
        filePattern = fullfile(scans_to_process_dir, 'raw_MRI-T1_index_*.csv');
        files = glob_sort_mtime(filePattern);
        raw_mris = readtable(files{1});

        % Create lists of raw processed MRIs
        idx = raw_mris.mri_processing_complete == 1;
        mri_dirs = raw_mris.mri_proc_dir(idx);
        mri_dirs = unique(mri_dirs);
        
        % Schedule PETs
        filePattern = fullfile(scans_to_process_dir, 'raw_PET_index_*.csv');
        files = glob_sort_mtime(filePattern);
        raw_pets = readtable(files{1});

        % Create the list of PET directories to process
        idx = raw_pets.pet_processing_complete == 1;
        pet_dirs = raw_pets.pet_proc_dir(idx);
        pet_dirs = unique(pet_dirs);
        
        run_dirs = [mri_dirs; pet_dirs];
     else
        run_dirs = abspath(cellvec(run_dirs));
        run_dirs = run_dirs(cellfun(@isfolder, run_dirs));
     end

    % Remove scans that have already been processed
    if ~overwrite
        run_dirs = run_dirs(cellfun(@(x) ~processed_wmap_exist(x), run_dirs));
    end 
    
    % Run W-map processing for scheduled files
    if isempty(run_dirs)
        fprintf('No scans to process\n');
        return;
    elseif length(run_dirs) == 1
        fprintf('Preparing to process 1 scan\n');
        cellfun(@(x) fprintf('  %s\n', x), run_dirs);
    else
        fprintf('Preparing to process %d scans\n', length(run_dirs));
        cellfun(@(x) fprintf('  %s\n', x), run_dirs);
    end

     % Process one scan
    if length(run_dirs) == 1
        process_single_wmap(run_dirs{1}, overwrite);
        
    % Process multiple scans in parallel
    else
        % Start a parallel pool
        n_workers = min(length(run_dirs), maxNumCompThreads);
        poolobj = parpool(n_workers);

        % Assign a worker to each scan
        parfor ii = 1:length(run_dirs)
            try
                process_single_wmap(run_dirs{ii}, overwrite);
            catch ME
                warning( ...
                    '\n\n\nERROR processing %s: %s\n\n\n', run_dirs{ii}, ...
                    getReport(ME, 'extended', 'hyperlinks', 'off') ...
                );
            end
        end

        % Close the parallel pool
        if ~isempty(poolobj)
            delete(poolobj);
        end
    end