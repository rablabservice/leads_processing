function process_mris(data_dir, log_dir, segment_brainstem, overwrite, verbose)
    % High-level function to select and process MRIs
    %
    % Parameters
    % ----------
    % data_dir:
    %   The directory that contains raw (unprocessed) MRI data in
    %   <data_dir>/raw and processed data in <data_dir>/processed
    % overwrite:
    %   If true, overwrite existing processed data
    % verbose:
    %   If true, print diagnostic information
    % ------------------------------------------------------------------
    arguments
        data_dir {mustBeFolder} = '/mnt/coredata/processing/leads/data'
        log_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/log'
        segment_brainstem logical = true
        overwrite logical = false
        verbose logical = true
    end

    % Format paths
    data_dir = abspath(data_dir);
    log_dir = abspath(log_dir);

    % Select MRIs to process
    mri_dirs = select_mris_to_process(log_dir, verbose);

    % Process MRIs in parallel
    parfor i = 1:length(mris_to_process)
        process_single_mri(mris_to_process(i), data_dir, segment_brainstem, overwrite, verbose);
    end
end
