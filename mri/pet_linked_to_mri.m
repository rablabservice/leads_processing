function linked_pet_dirs = pet_linked_to_mri(mri_dir)
    % Find all PET that have been linked to the MRI directory
    %
    % Parameters
    % ----------
    % mri_dir : char
    %     Path to the MRI directory of interest
    % ------------------------------------------------------------------
    arguments
        mri_dir char {mustBeNonempty}
    end

    linked_pet_dirs = {};

    % Get the parent directory of mri_dir
    parent_dir = fileparts(mri_dir);

    % List all subdirectories in the parent directory
    modalities = dir(parent_dir);

    % Filter out entries that are not directories and exclude MRI-T1
    modalities = modalities([modalities.isdir] & ~contains({modalities.name}, 'MRI-T1') ...
    & ~contains({modalities.name}, '.') & ~contains({modalities.name}, '..'));

    % Iterate through the modalities to check for an 'mri' symlink
    for i = 1:length(modalities)
        modality_path = fullfile(parent_dir, modalities(i).name);
        mri_symlink = fullfile(modality_path, 'mri');

        if isfolder(mri_symlink) 
            % Resolve the symlink target
            [status, symlink_target] = system(['readlink -f ', mri_symlink]);

            % Trim any trailing whitespace or newline characters
            symlink_target = abspath(strtrim(symlink_target));

            % Check if the symlink points to mri_dir
            if status == 0 && strcmp(symlink_target, abspath(mri_dir))
                linked_pet_dirs{end+1} = modality_path;
            end
        end
    end
end
