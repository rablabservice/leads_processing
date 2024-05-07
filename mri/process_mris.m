function process_mris(data_dir, log_dir, segment_brainstem, overwrite)
    % High-level function to select and process MRIs
    %
    % Parameters
    % ----------
    % data_dir:
    %   The directory that contains raw (unprocessed) MRI data in
    %   <data_dir>/raw and processed data in <data_dir>/processed
    % overwrite:
    %   If true, overwrite existing processed data
    % ------------------------------------------------------------------
    arguments
        data_dir {mustBeFolder} = '/mnt/coredata/processing/leads/data'
        log_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/log'
        segment_brainstem logical = true
        overwrite logical = false
    end

    % Format paths
    data_dir = abspath(data_dir);
    log_dir = abspath(log_dir);

    % Select MRIs to process
    mri_dirs = select_mris_to_process(log_dir);

    % Process MRIs in parallel
    parfor ii = 1:length(mri_dirs)
        process_single_mri(mri_dirs(ii), data_dir, segment_brainstem, overwrite);
    end
end
