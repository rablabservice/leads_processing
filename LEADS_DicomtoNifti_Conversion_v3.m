% Convert dicoms downloaded from LONI to nifti's 
% available in newdata folder and delete dicoms
% June 2023; nidhi.mundada@ucsf.edu
% Updated October 2023; piyush.maiti@ucsf.edu

% Converting MRIs 
fprintf(" * ------------ MRIs ------------")
cd (path_newmris)
checkAndConvertDICOM();

% Converting FBBs  
fprintf(" * ------------ FBB ------------")
cd (path_newfbbs)
checkAndConvertDICOM();

% Converting FTPs 
fprintf(" * ------------ FTP ------------")
cd (path_newftps)
checkAndConvertDICOM();

% Converting FDGs
fprintf(" * ------------ FDG ------------")
cd (path_newfdgs)
checkAndConvertDICOM();

% -----------------------------------------------      Function to Convert DICOM to NIfTi    -----------------------------------------------

function checkAndConvertDICOM()
    % This function checks for NIfTi and DICOM files in the current directory and its subdirectories.
    
    % Check for Existing NifTi and DICOM files
    nifiles = rdir('**/*.nii');
    dcmfiles = rdir('**/*.dcm');

    % If no NIfTi or DICOM files are found
    if isempty(dcmfiles) && isempty(nifiles)
        fprintf('\nNO DICOM or NIfTi Files available for %s\n', pwd);

    % If NIfTi files are found
    elseif ~isempty(nifiles)
        fprintf('\nNIfTi Files Found :\n');
        for i = 1:length(nifiles)
            fprintf('%s\n', nifiles(i).name);
        end
        fprintf('\nMoving on to the next Subject \n');

    % If DICOM files are found
    elseif ~isempty(dcmfiles)
        fprintf('DICOM Files found at %s', pwd);
        
        dcmFilePaths = {dcmfiles(:).name}';
        parentDirs = cellfun(@(x) fullfile(pwd, fileparts(x)), dcmFilePaths, 'UniformOutput', false);
        moddirs = unique(parentDirs);
        
        for i = 1:length(moddirs)
            cd(moddirs{i});

            fprintf('\nNow Converting: %s\n', moddirs{i});
            
            % dcm2niix path
            dcm2niixPath = '/mnt/coredata/Projects/LEADS/script_f7p1/service/dcm2niix';
            command = sprintf('"%s" -o "%s" "%s"', dcm2niixPath, pwd, pwd);

            % Run the system command
            [~, cmdOutput] = system(command);
            disp(cmdOutput);
            
            convdicoms = dir('*.nii');
            convdicomNames = {convdicoms(:).name}';

            % Rename nifti file according to the naming requirements 
            newnameParts = strsplit(moddirs{i}, filesep);
            newname = strcat(strjoin(newnameParts(end-3:end), '_'), '.nii');
            movefile(convdicomNames{1}, newname);

            delete('*.mat', '*.dcm');
            clear convdicoms dcmfiles newname i
        end
        fprintf(2,'%d Scans have been converted from dicom to nifti format. \n', length(moddirs));
        fprintf('\n');
    end
end