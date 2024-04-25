function varargout = get_freesurfer_files(mri_dir, fmt, segment_brainstem, fs_dir);
    % Locate and return paths to FreeSurfer files we want to work with
    %
    % The first argument out is a logical indicating whether all files
    % exist, and subsequent arguments out are paths to the files. This
    % function does not require mri_dir or any of the output files to
    % exist when called.
    %
    % Usage
    % -----
    % [~, nuf, aparcf, brainstemf] = get_freesurfer_files(mri_dir)
    % [all_exist, nu_mgzf, aparc_mgzf, brainstem_mgzf] = get_freesurfer_files(mri_dir, 'mgz', true)
    % [all_exist, nu_niif, aparc_niif, brainstem_niif] = get_freesurfer_files(mri_dir, 'nii', true)
    % [all_exist, nu_niif, aparc_niif] = get_freesurfer_files(mri_dir, 'nii', false)
    %
    % Parameters
    % ----------
    % mri_dir : char or str array
    %     Path to the processed MRI directory. Output files are written
    %     here. Unless fs_dir is passed, input files are assumed to be
    %     at <mri_dir>/freesurfer/mri/ when fmt is 'mgz'
    % fmt : char or str array, optional
    %     - If 'nii' (default), the FreeSurfer output files are expected
    %       to be in .nii format in <mri_dir>.
    %     - If 'mgz', the FreeSurfer output files are expected to be in
    %       .mgz format in <mri_dir>/freesurfer/mri.
    % segment_brainstem : logical, optional
    %     If true (default), check for and return the brainstem
    %     sublabels file
    % fs_dir : char or str array
    %     Path to the FreeSurfer directory where the output files from
    %     recon-all are located. If empty (default), it is assumed to be
    %     <mri_dir>/freesurfer.
    %
    % Returns
    % -------
    % all_exist : logical
    %     True if all output files exist, false otherwise
    % nu_niif : char or str array
    %     Path to the nu file
    % aparc_niif : char or str array
    %     Path to the aparc+aseg file
    % brainstem_niif : char or str array
    %     Path to the brainstem sublabels file
    % ------------------------------------------------------------------
    arguments
        mri_dir {mustBeText}
        fmt {mustBeMember(fmt, {'nii', 'mgz'})} = 'nii'
        segment_brainstem logical = true
        fs_dir {mustBeText} = ''
    end

    % Format paths
    mri_dir = abspath(mri_dir);
    if isempty(fs_dir)
        fs_dir = fullfile(mri_dir, 'freesurfer');
    end

    % Get paths to the nu and aparc+aseg files
    switch fmt
        case 'mgz'
            fs_mri_dir = fullfile(fs_dir, 'mri');
            nuf = fullfile(fs_mri_dir, 'nu.mgz');
            aparcf = fullfile(fs_mri_dir, 'aparc+aseg.mgz');
            outputs = {nuf, aparcf};
            if segment_brainstem
                brainstemf = fullfile(fs_mri_dir, 'brainstemSsLabels.v12.FSvoxelSpace.mgz');
                outputs{end+1} = brainstemf;
            end
        case 'nii'
            scan_tag = get_scan_tag(mri_dir);
            nuf = fullfile(mri_dir, strjoin({scan_tag, 'nu.nii'}, '_'));
            aparcf = fullfile(mri_dir, strjoin({scan_tag, 'aparc+aseg.nii'}, '_'));
            outputs = {nuf, aparcf};
            if segment_brainstem
                brainstemf = fullfile(mri_dir, strjoin({scan_tag, 'brainstem_sublabels.nii'}, '_'));
                outputs{end+1} = brainstemf;
            end
    end

    % Check if all files exist
    all_exist = all(cellfun(@(x) exist(x, 'file'), outputs));
    outputs = [{all_exist}, outputs];

    % Assign outputs
    [varargout{1:nargout}] = outputs{:};
end
