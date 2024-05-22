function pet_dirs = queue_pets_to_process(scans_to_process_dir)
    % Return a cell array of processed PET directories
    % for scans that are scheduled for processing
    % ------------------------------------------------------------------
    arguments
        scans_to_process_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/scans_to_process'
    end

    % Format paths
    scans_to_process_dir = abspath(scans_to_process_dir);

    % Print the welcome message
    fprintf('- Loading list of PET scans to process\n');

    % Load the most recent raw_pets table
    filePattern = fullfile(scans_to_process_dir, 'raw_PET_index_*.csv');
    files = glob_sort_mtime(filePattern);
    fprintf('  (%s)\n', files{1});
    raw_pets = readtable(files{1});

    % Create the list of PET directories to process
    idx = raw_pets.scheduled_for_processing == 1;
    pet_dirs = raw_pets.pet_proc_dir(idx);

    % Print summary statistics across tracers
    n_scans = height(raw_pets);
    n_subjs = numel(unique(raw_pets.subj));
    n_scans_processed = sum(raw_pets.pet_processing_complete);
    n_subjs_processed = numel(unique(raw_pets.subj(raw_pets.pet_processing_complete == 1)));
    n_scans_flagged = sum(raw_pets.flag);
    n_subjs_flagged = numel(unique(raw_pets.subj(raw_pets.flag == 1)));
    n_scans_to_process = sum(raw_pets.scheduled_for_processing);
    n_subjs_to_process = numel(unique(raw_pets.subj(raw_pets.scheduled_for_processing == 1)));
    n_scans_to_reprocess = sum(raw_pets.scheduled_for_processing(raw_pets.pet_processing_complete == 1));
    fprintf('  * %d PET scans from %d subjects in total\n', n_scans, n_subjs);
    fprintf('  * %d PET scans from %d subjects have been fully processed\n', n_scans_processed, n_subjs_processed);
    fprintf('  * %d PET scans from %d subjects are flagged with issues that preclude processing\n', n_scans_flagged, n_subjs_flagged);
    fprintf('  * %d PET scans from %d subjects are scheduled for processing\n', n_scans_to_process, n_subjs_to_process);
    if n_scans_to_reprocess > 0
        n_subjs_to_reprocess = numel(unique(raw_pets.subj((raw_pets.pet_processing_complete == 1) & (raw_pets.scheduled_for_processing == 1))));
        fprintf('    ...including %d already processed PET scans from %d subjects that will be reprocessed\n', n_scans_to_reprocess, n_subjs_to_reprocess);
    end
    fprintf('\n');

    % Print summary statistics for each tracer
    unique_tracers = unique(raw_pets.tracer);
    for i = 1:length(unique_tracers)
        tracer = unique_tracers{i};
        grp = raw_pets(strcmp(raw_pets.tracer, tracer), :);
        n_scans = height(grp);
        n_subjs = numel(unique(grp.subj));
        n_scans_processed = sum(grp.pet_processing_complete);
        n_subjs_processed = numel(unique(grp.subj(grp.pet_processing_complete == 1)));
        n_scans_flagged = sum(grp.flag);
        n_subjs_flagged = numel(unique(grp.subj(grp.flag == 1)));
        n_scans_to_process = sum(grp.scheduled_for_processing);
        n_subjs_to_process = numel(unique(grp.subj(grp.scheduled_for_processing == 1)));
        n_scans_to_reprocess = sum(grp.scheduled_for_processing(grp.pet_processing_complete == 1));
        fprintf('%s\n%s\n', tracer, repmat('-', 1, length(tracer)));
        fprintf('  * %d %s scans from %d subjects in total\n', n_scans, tracer, n_subjs);
        fprintf('  * %d %s scans from %d subjects have been fully processed\n', n_scans_processed, tracer, n_subjs_processed);
        fprintf('  * %d %s scans from %d subjects are flagged with issues that preclude processing\n', n_scans_flagged, tracer, n_subjs_flagged);
        fprintf('  * %d %s scans from %d subjects are scheduled for processing\n', n_scans_to_process, tracer, n_subjs_to_process);
        if n_scans_to_reprocess > 0
            n_subjs_to_reprocess = numel(unique(grp.subj((grp.pet_processing_complete == 1) & (grp.scheduled_for_processing == 1))));
            fprintf('    ...including %d already processed %s scans from %d subjects that will be reprocessed\n', n_scans_to_reprocess, tracer, n_subjs_to_reprocess);
        end
        fprintf('\n');
    end
end
