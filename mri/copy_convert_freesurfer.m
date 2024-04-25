function [nu_niif, aparc_niif, varargout] = copy_convert_freesurfer(mri_dir, segment_brainstem, overwrite, verbose)
    % Convert FreeSurfer files from .mgz to .nii and return nifti paths
    %
    % Nifti outputs are saved in the processed MRI directory.
    %
    % Usage
    % -----
    % [nuf, aparcf, brainstemf] = copy_convert_freesurfer(mri_dir)
    % [nuf, aparcf] = copy_convert_freesurfer(mri_dir, false)
    %
    % Parameters
    % ----------
    % mri_dir : char or str array
    %     Path to the processed MRI directory where converted files will be stored.
    % segment_brainstem : logical, optional
    %     If true (default), convert and return the brainstem .mgz file as .nii.
    % overwrite : logical, optional
    %     If true, overwrite existing .nii files (default is false).
    % verbose : logical, optional
    %     If true (default), print details of the conversion process to the console.
    %
    % Returns
    % -------
    % nu_niif : char or str array
    %     Path to the converted nu.nii file.
    % aparc_niif : char or str array
    %     Path to the converted aparc+aseg.nii file.
    % brainstem_niif : char or str array, optional
    %     Path to the converted brainstem.nii file, returned only if
    %     segment_brainstem is true.
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
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

    % Convert .mgz files to .nii
    if mgz_exist
        if overwrite || ~nii_exist
            if verbose
                fprintf('- Converting .mgz files to .nii\n')
            end

            % Convert the nu.mgz
            cmd_nu = append(mri_convert, nu_mgzf, ' ', nu_niif);
            if verbose
                fprintf('  $ %s\n', cmd_nu);
            end
            run_system_cmd(cmd_nu, true, false);

            % Convert the aparc+aseg.mgz
            cmd_aparc = append(mri_convert, aparc_mgzf, ' ', aparc_niif);
            if verbose
                fprintf('  $ %s\n', cmd_aparc);
            end
            run_system_cmd(cmd_aparc, true, false);

            % Convert the brainstem.mgz
            if segment_brainstem
                cmd_brainstem = append(mri_convert, brainstem_mgzf, ' ', brainstem_niif);
                if verbose
                    fprintf('  $ %s\n', cmd_brainstem);
                end
                run_system_cmd(cmd_brainstem, true, false);
            end
        end
    else
        fprintf('Missing one or more expected FreeSurfer outputs for %s', mri_dir);
    end

    % Format varargout
    if segment_brainstem
        varargout{1} = brainstem_niif;
    end
end
