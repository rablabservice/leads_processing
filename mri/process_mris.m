function process_mris(data_dir, overwrite, verbose)
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
        data_dir {mustBeText} = '/mnt/coredata/processing/leads/data'
        overwrite logical = false
        verbose logical = true
    end

    % Get the full, normalized path to the data directory
    data_dir = abspath(data_dir);

    % Select MRIs to process
    mris_to_process = select_mris_to_process(data_dir, overwrite, verbose);

    % Process each MRI
    parfor i = 1:length(mris_to_process)
        process_single_mri(mris_to_process(i), data_dir, overwrite, verbose);
    end
end
