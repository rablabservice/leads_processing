function copy_convert_freesurfer(mri_dir, segment_brainstem, overwrite, verbose)
    % Copy FreeSurfer files to the subject's processed mri directory
    % and convert them from .mgz to .nii
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeText}
        segment_brainstem logical = true
        overwrite logical = false
        verbose logical = true
    end

    % Format inputs
    mri_dir = abspath(mri_dir);

    % Define the mri_convert command with appropriate flags
    mri_convert = 'mri_convert -it mgz -ot nii --out_orientation RAS ';

    % Get the scan tag and path to the freesurfer mri directory
    if segment_brainstem
        [mgz_exist, nu_mgzf, aparc_mgzf, brainstem_mgzf] = get_freesurfer_files(mri_dir, 'mgz', segment_brainstem);
        [nii_exist, nu_niif, aparc_niif, brainstem_niif] = get_freesurfer_files(mri_dir, 'nii', segment_brainstem);
    else
        [mgz_exist, nu_mgzf, aparc_mgzf] = get_freesurfer_files(mri_dir, 'mgz', segment_brainstem);
        [nii_exist, nu_niif, aparc_niif] = get_freesurfer_files(mri_dir, 'nii', segment_brainstem);
    end


    if mgz_exist
        if overwrite || ~nii_exist
            if verbose
                fprintf('Converting FreeSurfer files from .mgz to .nii for %s\n', mri_dir);
            end
            system(char(append(mri_convert, nu_mgzf, ' ', nu_niif)));
            system(char(append(mri_convert, aparc_mgzf, ' ', aparc_niif)));
            if segment_brainstem
                system(char(append(mri_convert, brainstem_mgzf, ' ', brainstem_niif)));
            end
        end
    else
        error('Missing one or more expected FreeSurfer outputs for %s', mri_dir);
    end
end
