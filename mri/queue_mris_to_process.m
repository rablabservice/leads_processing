function mri_dirs = queue_mris_to_process(scans_to_process_dir)
    % Return cell arrays of raw MRI files and processed MRI directories
    % for scans that are scheduled for processing
    % ------------------------------------------------------------------
    arguments
        scans_to_process_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/scans_to_process'
    end

    % Format paths
    scans_to_process_dir = abspath(scans_to_process_dir);

    % Print the welcome message
    fprintf('- Loading list of MRIs to process\n');

    % Load the most recent raw_MRIs table
    filePattern = fullfile(scans_to_process_dir, 'raw_MRI_index_*.csv');
    files = glob_sort_mtime(filePattern);
    fprintf('  (%s)\n', files{1});
    raw_mris = readtable(files{1});

    % Create lists of raw MRI files to process and their corresponding
    % processed MRI directories
    idx = raw_mris.scheduled_for_processing == 1;
    mri_dirs = raw_mris.mri_proc_dir(idx);

    % Print summary statistics
    n_scans = height(raw_mris);
    n_subjs = numel(unique(raw_mris.subj));
    n_scans_freesurfer_complete = sum(raw_mris.freesurfer_complete);
    n_subjs_freesurfer_complete = numel(unique(raw_mris.subj(raw_mris.freesurfer_complete == 1)));
    n_scans_processed = sum(raw_mris.mri_processing_complete);
    n_subjs_processed = numel(unique(raw_mris.subj(raw_mris.mri_processing_complete == 1)));
    n_scans_to_process = sum(raw_mris.scheduled_for_processing);
    n_subjs_to_process = numel(unique(raw_mris.subj(raw_mris.scheduled_for_processing == 1)));
    n_baseline_scans_to_process = sum(raw_mris.scheduled_for_processing(raw_mris.mri_scan_number == 0));
    n_followup_scans_to_process = sum(raw_mris.scheduled_for_processing(raw_mris.mri_scan_number > 0));
    n_scans_to_reprocess = sum(raw_mris.scheduled_for_processing(raw_mris.mri_processing_complete == 1));
    fprintf('  * %d MRIs from %d subjects in total\n', n_scans, n_subjs);
    fprintf('  * %d MRIs from %d subjects have completed FreeSurfer processing\n', n_scans_freesurfer_complete, n_subjs_freesurfer_complete);
    fprintf('  * %d MRIs from %d subjects have been fully processed\n', n_scans_processed, n_subjs_processed);
    fprintf('  * %d MRIs from %d subjects are scheduled for processing\n', n_scans_to_process, n_subjs_to_process);
    fprintf('    - including %d baseline MRIs\n', n_baseline_scans_to_process);
    fprintf('    - and       %d follow-up MRIs\n', n_followup_scans_to_process);
    if n_scans_to_reprocess > 0
        n_subjs_to_reprocess = numel(unique(raw_mris.subj((raw_mris.mri_processing_complete == 1) & (raw_mris.scheduled_for_processing == 1))));
        fprintf('    ...including %d already processed MRIs from %d subjects that will be reprocessed\n', n_scans_to_reprocess, n_subjs_to_reprocess);
    end
    fprintf('\n');

    % Print the full list of scans to process, separating each scan by a
    % space and adding a newline when adding the next scan causes
    % the current line length to exceed max_line
    max_line = 115;
    current_line = '';
    scan_tags = cellfun(@get_scan_tag, mri_dirs, 'UniformOutput', false);
    fprintf('All MRIs scheduled to be processed:\n');
    for i = 1:length(scan_tags)
        new_entry = scan_tags{i};
        if isempty(current_line)
            new_line = new_entry;
        else
            new_line = append(current_line, '  ', new_entry);
        end

        % Check if the new line length exceeds max_line
        if length(new_line) > max_line
            % Print the current line and start a new one
            fprintf('%s\n', current_line);
            current_line = new_entry;
        else
            % Update the current line
            current_line = new_line;
        end
    end

    % Print any remaining text in the current line
    if ~isempty(current_line)
        fprintf('%s\n', current_line);
    end
    fprintf('\n');
end
