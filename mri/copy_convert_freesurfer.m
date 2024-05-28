function outfiles = copy_convert_freesurfer(mri_dir, fid, overwrite)
    % Convert FreeSurfer files from .mgz to .nii and return nifti paths
    %
    % Converted .nii files are saved in mri_dir
    %
    % Parameters
    % ----------
    % mri_dir : char or str array
    %     Path to the processed MRI directory where converted files will be stored
    % fid : int, optional
    %     File identifier for logging (default is 1 for stdout)
    % overwrite : logical, optional
    %     If true, overwrite existing .nii files (default is false)
    %
    % Returns
    % -------
    % outfiles : struct
    %     Struct array with paths to each converted nifti file (if it
    %     exists):
    %     - nu    : The nu file
    %     - aparc : The aparc+aseg file
    %     - bstem : The brainstem sublabels file
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
        fid {mustBeNumeric} = 1
        overwrite logical = false
    end

    % Format inputs
    mri_dir = abspath(mri_dir);

    % Define the mri_convert command with appropriate flags
    mri_convert = 'mri_convert -it mgz -ot nii --out_orientation RAS';

    % Get the scan tag and path to the freesurfer mri directory
    fs_mgzfs = get_freesurfer_files(mri_dir, 'mgz');
    fs_niifs = get_freesurfer_files(mri_dir, 'nii');

    % Loop over each file in fs_mgzfs, check if the corresponding file
    % in fs_niftis exists, and convert .mgz to .nii if necessary
    log_append(fid, '- Converting FreeSurfer .mgz files to .nii');
    fs_fields = fieldnames(fs_mgzfs);
    for ii = 1:length(fs_fields)
        mgzf = fs_mgzfs.(fs_fields{ii});
        niif = fs_niifs.(fs_fields{ii});
        if isfile(mgzf)
            if overwrite || ~isfile(niif)
                cmd = sprintf('%s %s %s', mri_convert, mgzf, niif);
                msg = sprintf( ...
                    '  * %s ->\n              %s', ...
                    basename(mgzf), ...
                    basename(niif) ...
                );
                log_append(fid, msg);
                run_system(cmd, 1, false, true);
            else
                msg = sprintf('  * %s exists, will not overwrite', basename(niif));
                log_append(fid, msg);
            end
        end
        % Keep only outfiles that exist
        if isfile(niif)
            outfiles.(fs_fields{ii}) = niif;
        end
    end
end
