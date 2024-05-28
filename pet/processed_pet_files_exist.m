function all_exist = processed_pet_files_exist(pet_dir)
    % Return true if all freesurfer directory is complete
    %
    % Parameters
    % ----------
    % pet_dir : char or str
    %     The directory that contains the processed MRI data
    %
    % Returns
    % -------
    % all_exist : logical
    %    true if all processed MRI files are present, else false
    % ------------------------------------------------------------------
    arguments
        pet_dir {mustBeText} = ''
    end

    pet_files = get_processed_pet_files(pet_dir);

    all_exist = all(isfile(cellvec(pet_files)));
