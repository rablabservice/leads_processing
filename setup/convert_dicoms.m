function convert_dicoms(topdir)
    % Convert .dcm files to .nii for all nested directories in
    % topdir with dicoms but no existing nifti.
    % ------------------------------------------------------------------
    arguments
        topdir {mustBeFolder}
    end

    % Format paths
    topdir = abspath(topdir);

    % Get path to dcm2niix
    dcm2niix = '/home/mac/dschonhaut/bin/dcm2niix -d 9';

    % Find all dicoms in raw
    dcm_files = dir(fullfile(topdir, '**', '*.dcm'));
    dcm_paths = fullfile({dcm_files.folder}', {dcm_files.name}');
    dcm_dirs = unique(cellfun(@fileparts, dcm_paths, 'UniformOutput', false));

    % Find all niftis in raw
    nii_files = dir(fullfile(topdir, '**', '*.nii*'));
    nii_paths = fullfile({nii_files.folder}', {nii_files.name}');
    nii_dirs = unique(cellfun(@fileparts, nii_paths, 'UniformOutput', false));

    % Find directories with dicoms but no nifti
    conv_dirs = setdiff(dcm_dirs, nii_dirs);
    n_conv = length(conv_dirs);
    fprintf('  * Found %d scan directories with dicoms to convert\n', n_conv);

    % Convert dicoms to nifti, parallelizing if there are more than 10
    % directories with dicoms to convert
    if n_conv < 10
        for i = 1:n_conv
            scan_dir = abspath(conv_dirs{i});
            cmd = char(append(dcm2niix, ' -o ', scan_dir, ' ', scan_dir));
            run_system(cmd, 1, false, true);
        end
    else
        parfor i = 1:n_conv
            scan_dir = abspath(conv_dirs{i});
            cmd = char(append(dcm2niix, ' -o ', scan_dir, ' ', scan_dir));
            run_system(cmd, 1, false, true);
        end
    end
end
