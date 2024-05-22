function all_exist = freesurfer_files_exist(mri_dir, fs_dir, segment_brainstem)
    % Return true if all expected freesurfer *.mgz files exist
    %
    % Parameters
    % ----------
    % mri_dir : char|string
    %     Path to the directory where FreeSurfer converted .nii files
    %     are located. If fs_dir is provided, mri_dir is ignored
    % fs_dir : char|string
    %     Path to the FreeSurfer directory created by recon-all. If this
    %     field is empty, it is assumed to be <mri_dir>/freesurfer
    %
    % Returns
    % -------
    % all_exist : logical
    %    true if all FreeSurfer files are present, else false. The
    %    following files are checked:
    %    - nu.mgz
    %    - aparc+aseg.mgz
    %    - brainstemSsLabels.v12.FSvoxelSpace.mgz
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeText} = ''
        fs_dir {mustBeText} = ''
        segment_brainstem logical = true
    end

    % Get paths to the FreeSurfer .mgz files
    fsfs = get_freesurfer_files(mri_dir, 'mgz', fs_dir);
    if ~segment_brainstem
        rmfield(fsfs, 'bstem');
    end

    % Check if all files exist
    all_exist = all(isfile(cellvec(fsfs)));
