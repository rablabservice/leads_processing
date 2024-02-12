function output = getScanInfo(subjDir)
% Return info on each scan in subject directory.
subj = fileparts(subjDir);
dicoms = dir(fullfile(subjDir, '**', '*.dcm'));
niftis = dir(fullfile(subjDir, '**', '*.nii*'));
scanDirs = unique({dicoms.folder, niftis.folder});
output = [];

for i = 1:length(scanDirs)
    scanDir = scanDirs{i};
    petDate = getAcqDate(scanDir);
    notes = "";
    if ~isempty(petDate)
        niftis = dir(fullfile(scanDir, '*.nii*'));
        dicoms = dir(fullfile(scanDir, '*.dcm'));
        if isempty(niftis) && ~isempty(dicoms)
            niftis = dcm2niix(scanDir); % Assuming dcm2niix is a custom or installed function
        end
        if length(niftis) == 1
            niiPath = fullfile(niftis(1).folder, niftis(1).name);
        elseif isempty(niftis)
            niiPath = "";
            notes = notes + "Can't find reconstructed nifti for " + scanDir + "; ";
        else
            niiPath = "";
            notes = notes + "Multiple niftis found for " + scanDir + "; ";
        end
        if ~isempty(niiPath)
            tracer = getTracer(niiPath);
            if isempty(tracer)
                notes = notes + "Can't parse PET tracer from " + niiPath + "; ";
            end
            inputRes = getInputRes(niiPath);
            if isempty(inputRes)
                notes = notes + "Can't parse starting resolution from " + niiPath + "; ";
            end
        else
            tracer = "";
            inputRes = [];
        end
    else
        niiPath = "";
        tracer = "";
        inputRes = [];
        notes = notes + "Can't find PET acquisition date for " + scanDir + "; ";
    end
    output = [output; {subj, petDate, tracer, inputRes, niiPath, notes}];
end

% Convert output to table for easier handling, similar to pandas DataFrame
output = cell2table(output, 'VariableNames', {"subj", "petDate", "tracer", "inputRes", "rawPetf", "notes"});
end

function acqDate = getAcqDate(filepath)
% Return the acquisition date as YYYY-MM-DD.
parts = strsplit(filepath, filesep);
for i = length(parts):-1:1
    d = parts{i};
    try
        datetime(d(1:10), 'InputFormat', 'yyyy-MM-dd');
        acqDate = d(1:10);
        return;
    catch
    end
end
acqDate = [];
end

function tracer = getTracer(filepath)
% Return the PET tracer used from filepath to the recon'd nifti.
[~, basename, ~] = fileparts(filepath);
basename = lower(basename);
tracers = {'fbb', 'florbetaben', 'neuraceq', 'fbp', 'florbetapir', 'av45', 'av-45', 'amyvid', 'flutafuranol', 'nav4694', 'nav-4694', 'azd4694', 'azd-4694', 'pib', 'pittsburgh compound b', 'pittsburgh compound-b'};
tracerNames = {'FBB', 'FBB', 'FBB', 'FBP', 'FBP', 'FBP', 'FBP', 'FBP', 'NAV', 'NAV', 'NAV', 'NAV', 'NAV', 'PIB', 'PIB', 'PIB'};
for i = 1:length(tracers)
    if contains(basename, tracers{i})
        tracer = tracerNames{i};
        return;
    end
end
tracer = [];
end

function inputRes = getInputRes(filepath)
% Get the starting resolution from file basename.
[~, basename, ~] = fileparts(filepath);
basename = lower(basename);
if contains(basename, 'uniform_6mm_res')
    inputRes = [6, 6, 6];
elseif contains(basename, 'uniform_8mm_res')
    inputRes = [8, 8, 8];
elseif contains(basename, 'uniform_') && contains(basename, 'mm_res')
    resStr = extractBetween(basename, 'uniform_', 'mm_res');
    resNum = str2double(resStr);
    inputRes = [resNum, resNum, resNum];
else
    inputRes = [];
end
end
