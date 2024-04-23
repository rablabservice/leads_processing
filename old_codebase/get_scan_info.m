function output = get_scan_info(subjDir, scanTypeMapFile)
    % Parse subject dir and return subject ID, scan date, and scan type.

    % Load the file with modality input to output name mappings
    switch nargin
        case 1
            scanTypeMapFile = '/mnt/coredata/processing/leads/metadata/ssheets/scan_types_and_tracers.csv';
    end
    scanTypeMap = readtable(scanTypeMapFile, 'TextType', 'string');

    % Find all unique nifti-containing subdirectories
    scanDirs = unique({dir(fullfile(subjDir, '**', '*.nii*')).folder});
    if isempty(scanDirs)
        output = cell2table(cell(0, 4), 'VariableNames', {'subj', 'scanDate', 'scanType', 'scanPath'});
        return;
    end

    % Get the subject ID from the base directory name
    subj = getSubjectID(subjDir)

    output = [];
    for i = 1:length(scanDirs)
        scanDir = scanDirs{i};
        scanDate = getScanDate(scanDir);
        niftis = dir(fullfile(scanDir, '*.nii*'));
        niiPath = fullfile(niftis(1).folder, niftis(1).name);
        scanType = getScanType(niiPath, scanTypeMap);
        output = [output; {subj, scanDate, scanType, niiPath}];
    end

    % Convert output to table for easier handling, similar to pandas DataFrame
    output = cell2table(output, 'VariableNames', {'subj', 'scanDate', 'scanType', 'scanPath'});
end

function subj = getSubjectID(subjDir)
    % Return the subject ID from the base directory name.
    [~, subj, ~] = fileparts(subjDir);
end

function scanDate = getScanDate(filepath)
    % Return the scan acquisition date as YYYY-MM-DD.
    parts = strsplit(filepath, filesep);
    for i = length(parts):-1:1
        d = parts{i};
        try
            datetime(d(1:10), 'InputFormat', 'yyyy-MM-dd');
            scanDate = d(1:10);
            return;
        catch
        end
    end
    scanDate = '';
end

function scanType = getScanType(filepath, scanTypeMap)
    % Parse the filepath and return the scan type.
    [~, basename, ~] = fileparts(filepath);
    basename = lower(string(basename));

    % Search through each scan type in the mapping
    for i = 1:height(scanTypeMap)
        if contains(basename, scanTypeMap.name_in{i})
            scanType = scanTypeMap.name_out{i};
            return;
        end
    end

    % If no scan type is found, return an empty string
    scanType = '';
end
