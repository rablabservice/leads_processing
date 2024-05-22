function all_exist = processed_mri_files_exist(mri_dir)
    % Return true if all freesurfer directory is complete
    %
    % Parameters
    % ----------
    % mri_dir : char or str
    %     The directory that contains the processed MRI data
    %
    % Returns
    % -------
    % all_exist : logical
    %    true if all processed MRI files are present, else false
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeText} = ''
    end

    mrifs = get_processed_mri_files(mri_dir);

    all_exist = all(isfile(cellvec(mrifs)));
