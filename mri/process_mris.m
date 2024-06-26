function process_mris( ...
    mri_dirs, ...
    scans_to_process_dir, ...
    overwrite, ...
    segment_brainstem, ...
    process_freesurfer, ...
    process_post_freesurfer ...
)
    % Process all MRIs that are scheduled for processing
    % in the latest raw_MRI_index file
    %
    % Parameters
    % ----------
    % mri_dirs : cell
    %     List of MRI directories to process. This variable overrides
    %     the default behavior of selecting scans to process based on
    %     the latest raw_MRI_index file in scans_to_process_dir
    % scans_to_process_dir : char or str
    %     The directory that stores raw_MRI_index files
    % overwrite : logical
    %     If true, overwrite existing processed data
    % segment_brainstem : logical
    %     If true, segment the brainstem using segmentBS.sh
    % process_freesurfer : logical
    %     If true, run recon-all on the raw MRI
    % process_post_freesurfer : logical
    %     If true, run post-FreeSurfer processing
    % ------------------------------------------------------------------
    arguments
        mri_dirs = {}
        scans_to_process_dir {mustBeText} = '/mnt/coredata/processing/leads/metadata/scans_to_process'
        overwrite logical = false
        segment_brainstem logical = true
        process_freesurfer logical = true
        process_post_freesurfer logical = true
    end

    % Format paths
    scans_to_process_dir = abspath(scans_to_process_dir);

    % Load the list of MRIs to process
    fprintf('Checking list of submitted scans to process...\n');
    if isempty(mri_dirs)
        mri_dirs = queue_mris_to_process(scans_to_process_dir);
    else
        mri_dirs = abspath(cellvec(mri_dirs));
        mri_dirs = mri_dirs(cellfun(@isfolder, mri_dirs));
    end
    raw_mrif = '';

    % Remove scans that have already been processed
    if ~overwrite
        mri_dirs = mri_dirs(cellfun(@(x) ~processed_pet_files_exist(x), mri_dirs));
    end

    % Are there are any scans to process?
    if isempty(mri_dirs)
        fprintf('No MRIs to process\n');
        return;
    elseif length(mri_dirs) == 1
        fprintf('Preparing to process 1 MRI\n');
        cellfun(@(x) fprintf('  %s\n', x), mri_dirs);
    else
        fprintf('Preparing to process %d MRIs\n', length(mri_dirs));
        cellfun(@(x) fprintf('  %s\n', x), mri_dirs);
    end

    % Process one scan
    if length(mri_dirs) == 1
        process_single_mri( ...
            mri_dirs{1}, ...
            overwrite, ...
            raw_mrif, ...
            segment_brainstem, ...
            process_freesurfer, ...
            process_post_freesurfer ...
        );
    % Process multiple scans in parallel
    else
        % Start a parallel pool
        n_workers = min(length(mri_dirs), maxNumCompThreads);
        poolobj = parpool(n_workers);

        % Assign a worker to each scan
        parfor ii = 1:length(mri_dirs)
            try
                process_single_mri( ...
                    mri_dirs{ii}, ...
                    overwrite, ...
                    raw_mrif, ...
                    segment_brainstem, ...
                    process_freesurfer, ...
                    process_post_freesurfer ...
                );
            catch ME
                warning( ...
                    '\n\n\nERROR processing %s: %s\n\n\n', mri_dirs{ii}, ...
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
