function process_mris( ...
    overwrite, ...
    scans_to_process_dir, ...
    segment_brainstem, ...
    process_freesurfer, ...
    process_post_freesurfer ...
)
    % Process all MRIs that are scheduled for processing
    % in the latest log file
    %
    % Parameters
    % ----------
    % overwrite : logical
    %     If true, overwrite existing processed data
    % scans_to_process_dir : char or str
    %     The directory that stores log files
    % segment_brainstem : logical
    %     If true, segment the brainstem using segmentBS.sh
    % process_freesurfer : logical
    %     If true, run recon-all on the raw MRI
    % process_post_freesurfer : logical
    %     If true, run post-FreeSurfer processing
    % ------------------------------------------------------------------
    arguments
        overwrite logical = false
        scans_to_process_dir {mustBeFolder} = '/mnt/coredata/processing/leads/metadata/scans_to_process'
        segment_brainstem logical = true
        process_freesurfer logical = true
        process_post_freesurfer logical = true
    end

    % Format paths
    scans_to_process_dir = abspath(scans_to_process_dir);

    % Load the list of MRIs to process
    mri_dirs = queue_mris_to_process(scans_to_process_dir);
    raw_mrif = '';

    % Process MRIs in parallel
    parfor ii = 1:length(mri_dirs)
        try
            process_single_mri( ...
                mri_dirs{ii}, ...
                raw_mrif, ...
                overwrite, ...
                segment_brainstem, ...
                process_freesurfer, ...
                process_post_freesurfer ...
            );
        catch ME
            warning('\n\n\nERROR processing %s: %s\n\n\n', mri_dirs{ii}, ME.message);
        end
    end
end
