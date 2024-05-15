function process_mris(overwrite, log_dir, segment_brainstem)
    % High-level function to select and process MRIs
    %
    % Parameters
    % ----------
    % overwrite : logical
    %     If true, overwrite existing processed data
    % data_dir : char or str
    %     The directory that contains raw (unprocessed) MRI data in
    %     <data_dir>/raw and processed data in <data_dir>/processed
    % log_dir : char or str
    %     The directory that stores log files
    % segment_brainstem : logical
    %     If true, segment the brainstem using segmentBS.sh
    % ------------------------------------------------------------------
    arguments
        overwrite logical = false
        log_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/log'
        segment_brainstem logical = true

    end

    % Format paths
    log_dir = abspath(log_dir);

    % Select MRIs to process
    [raw_mrifs, mri_dirs] = select_mris_to_process(log_dir);

    % Process MRIs in parallel
    parfor ii = 1:length(mri_dirs)
        process_single_mri(raw_mrifs(ii), mri_dirs(ii), overwrite, segment_brainstem);
    end
end
