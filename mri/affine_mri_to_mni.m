function affine_mri_to_mni(mri_dir, overwrite, verbose)
    % Calculate the full affine transform from native MRI to MNI space
    %
    % Parameters
    % ----------
    % mri_dir : char or str array
    %   The directory that contains the processed MRI data
    % overwrite : logical, optional
    %   If true, overwrite existing files
    % verbose : logical, optional
    %   If true, print diagnostic information
    %
    % Files created
    % -------------
    % - <mri_dir>/a<scan_tag>_nu.nii
    % - <mri_dir>/<scan_tag>_affine-to-mni.mat (or some file like this that we can use
    %                                           later on PET files coreg'd to native MRI)
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeFolder}
        overwrite logical = false
        verbose logical = true
    end

    % Format paths
    mri_dir = abspath(mri_dir);
    scan_tag = get_scan_tag(mri_dir);
    [~, nuf] = get_freesurfer_files(mri_dir, 'nii');

    % Check if affine transform has already been calculated
    wnuf = fullfile(mri_dir, append('a', basename(nuf)));
    if exist(anuf, 'file') && ~overwrite
        return
    end

    % Calculate affine transform to MNI space
    if verbose
        fprintf('- Segmenting nu.nii\n');
    end
    % DO STUFF
end
