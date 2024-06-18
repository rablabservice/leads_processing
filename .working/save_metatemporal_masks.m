% Define path to the processed data directory
proc_dir = '/mnt/coredata/processing/leads/data/processed';

% Get a list of all directories one level down from proc_dir
subj_dirs = dir(fullfile(proc_dir, '*'));
subj_dirs = subj_dirs([subj_dirs.isdir]); % Keep only directories

% Initialize an empty cell array to store the desired directories
mri_dirs = {};

% Loop through each directory in level 1
for ii = 1:length(subj_dirs)
    if strcmp(subj_dirs(ii).name, '.') || strcmp(subj_dirs(ii).name, '..')
        continue;
    end

    % Get the full path of the level 1 directory
    level1_path = fullfile(proc_dir, subj_dirs(ii).name);

    % Get a list of all directories one level down from the current level 1 directory
    scan_dirs = dir(fullfile(level1_path, '*'));
    scan_dirs = scan_dirs([scan_dirs.isdir]); % Keep only directories

    % Loop through each directory in level 2
    for jj = 1:length(scan_dirs)
        if strcmp(scan_dirs(jj).name, '.') || strcmp(scan_dirs(jj).name, '..')
            continue;
        end

        % Check if the directory name starts with 'MRI-T1'
        if startsWith(scan_dirs(jj).name, 'MRI-T1')
            % Add the directory path to the result cell array
            mri_dirs{end+1} = fullfile(level1_path, scan_dirs(jj).name); %#ok<SAGROW>
        end
    end
end

% Convert the result cell array to a cellstr array
mri_dirs = cellvec(mri_dirs);

% Get paths to the aparcs
aparc_files = cellfun(@(x) get_freesurfer_files(x).aparc, mri_dirs, 'UniformOutput', false);

% Display the result
fprintf('Found %d aparc+aseg.nii files\n', length(aparc_files));

% Start a parallel pool
max_workers = 16;
n_workers = min(length(aparc_files), max_workers);
poolobj = parpool(n_workers);

% Assign a worker to each scan
parfor ii = 1:length(aparc_files)
    try
        save_mask_metatemporal(aparc_files{ii});
    catch ME
        warning( ...
            '\n\n\nERROR processing %s: %s\n\n\n', aparc_files{ii}, ...
            getReport(ME, 'extended', 'hyperlinks', 'off') ...
        );
    end
end

% Close the parallel pool
if ~isempty(poolobj)
    delete(poolobj);
end
