function convert_dicoms(newdata_dir)
    % Convert .dcm files to .nii for all nested directories in newdata_dir
    % with dicoms but no existing nifti.
    % ------------------------------------------------------------------
    arguments
        newdata_dir {mustBeFolder}
    end

    % Format paths
    newdata_dir = abspath(newdata_dir);

    % If newdata_dir is empty, there's nothing to do here
    if isempty(dir(newdata_dir))
        fprintf('- No newdata to convert to nifti\n');
        return
    end

    % Get path to dcm2niix
    dcm2niix = '/home/mac/dschonhaut/bin/dcm2niix -d 9';

    % Print the welcome message
    fprintf('\n- Converting newdata dicoms to nifti\n');

    % Find all dicoms in raw
    dcm_files = dir(fullfile(newdata_dir, '**', '*.dcm'));
    dcm_paths = fullfile({dcm_files.folder}', {dcm_files.name}');
    dcm_dirs = unique(cellfun(@fileparts, dcm_paths, 'UniformOutput', false));

    % Find all niftis in raw
    nii_files = dir(fullfile(newdata_dir, '**', '*.nii*'));
    nii_paths = fullfile({nii_files.folder}', {nii_files.name}');
    nii_dirs = unique(cellfun(@fileparts, nii_paths, 'UniformOutput', false));

    % Find directories with dicoms but no nifti
    conv_dirs = setdiff(dcm_dirs, nii_dirs);
    n_conv = length(conv_dirs);
    fprintf('  * %d directories have dicoms but no nifti\n', n_conv);

    % Convert dicoms to nifti
    for i = 1:n_conv
        scan_dir = abspath(conv_dirs{i});
        cmd = char(append(dcm2niix, ' -o ', scan_dir, ' ', scan_dir));
        run_system_cmd(cmd);
    end
end
