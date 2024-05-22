function outfiles = get_freesurfer_files(mri_dir, fmt, fs_dir);
    % Locate and return paths to FreeSurfer files we want to work with
    %
    % This function does not require mri_dir or any of the output files
    % to exist when called.
    %
    % Usage
    % -----
    % fs_niifs = get_freesurfer_files(mri_dir)
    % fs_mgzfs = get_freesurfer_files(mri_dir, 'mgz')
    %
    % Parameters
    % ----------
    % mri_dir : char|string
    %     Path to the directory where FreeSurfer converted .nii files
    %     are located
    % fmt : char|string, optional
    %     - If 'nii' (default), this function returns paths to the
    %       FreeSurfer converted .nii files in <mri_dir>
    %     - If 'mgz', the FreeSurfer output files are expected to be in
    %       .mgz format in <mri_dir>/freesurfer/mri.
    % fs_dir : char|string
    %     Path to the FreeSurfer directory created by recon-all. If this
    %     field is empty (default), it is assumed to be
    %     <mri_dir>/freesurfer
    %
    % Returns
    % -------
    % outfiles : struct
    %     Struct array with paths to:
    %     - nu    : The nu file
    %     - aparc : The aparc+aseg file
    %     - bstem : The brainstem sublabels file
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeText} = ''
        fmt {mustBeMember(fmt, {'nii', 'mgz'})} = 'nii'
        fs_dir {mustBeText} = ''
    end

    % Format parameters
    mri_dir = abspath(mri_dir);

    % Get paths to the nu and aparc+aseg files
    switch fmt
        case 'mgz'
            if isempty(fs_dir)
                fs_dir = fullfile(mri_dir, 'freesurfer');
            end
            fs_mri_dir = fullfile(fs_dir, 'mri');
            outfiles.nu = fullfile(fs_mri_dir, 'nu.mgz');
            outfiles.aparc = fullfile(fs_mri_dir, 'aparc+aseg.mgz');
            outfiles.bstem = fullfile(fs_mri_dir, 'brainstemSsLabels.v12.FSvoxelSpace.mgz');
        case 'nii'
            scan_tag = get_scan_tag(mri_dir);
            outfiles.nu = fullfile(mri_dir, append(scan_tag, '_nu.nii'));
            outfiles.aparc = fullfile(mri_dir, append(scan_tag, '_aparc+aseg.nii'));
            outfiles.bstem = fullfile(mri_dir, append(scan_tag, '_brainstem_sublabels.nii'));
    end
end
