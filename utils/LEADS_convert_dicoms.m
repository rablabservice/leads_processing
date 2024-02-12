%% Recursively search the newdata directory, find directories containing
%% DICOM files but no NIfTI files, and convert DICOMs to NIfTI in these
%% directories
dcm2niix = fullfile(dirs('code'), 'utils', 'dcm2niix');

% Find all dicoms in newdata
dcm_files = dir(fullfile(dirs('newdata'), '**', '*.dcm'));
dcm_paths = fullfile({dcm_files.folder}', {dcm_files.name}');
dcm_dirs = unique(cellfun(@fileparts, dcm_paths, 'UniformOutput', false));

% Find all niftis in newdata
nii_files = dir(fullfile(dirs('newdata'), '**', '*.nii*'));
nii_paths = fullfile({nii_files.folder}', {nii_files.name}');
nii_dirs = unique(cellfun(@fileparts, nii_paths, 'UniformOutput', false));

% Find directories with dicoms but no niftis
conv_dirs = setdiff(dcm_dirs, nii_dirs);

% Convert dicoms to niftis
for i = 1:length(conv_dirs)
    cd(conv_dirs{i});
    cmd = sprintf('"%s" -o "%s" "%s"', dcm2niix, conv_dirs{i}, conv_dirs{i});
    [status, cmdout] = system(cmd);
end
