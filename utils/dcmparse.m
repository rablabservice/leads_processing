function dicom_fields = dcmparse(dicom_dir, dicom_fields)
    % Return a table that maps each dicom field to an array of unique values across dicom headers in dicom_dir files
    arguments
        dicom_dir string = pwd()
        dicom_fields cell = { ...
            'SeriesDescription', 'SeriesDate', 'SeriesTime', 'NumberOfSlices', 'NumberOfTimeSlices', 'AcquisitionTime', 'ContentTime', 'FrameReferenceTime', 'ActualFrameDuration', 'SeriesNumber', 'AcquisitionNumber', 'InstanceNumber', 'ImageIndex', 'SliceLocation', 'AxialAcceptance', 'DecayFactor',  'Filename' ...
        }
    end

    % Find all DICOM files in the specified directory
    alldicoms = dir(fullfile(dicom_dir, "*.IMA"));

    % Initialize a struct to store unique field values
    unique_fields = struct();

    % Loop through each DICOM file
    for i = 1:length(alldicoms)
        file_path = fullfile(alldicoms(i).folder, alldicoms(i).name);

        % Try to read DICOM info, skip if failure occurs (file may not be a valid DICOM)
        try
            info = dicominfo(file_path);
        catch
            continue;
        end

        % Loop through each requested field
        for field = dicom_fields
            if isfield(info, field{1})
                % Store unique values of each field in a cell array
                if isfield(unique_fields, field{1})
                    unique_fields.(field{1}) = [unique_fields.(field{1}), info.(field{1})];
                else
                    unique_fields.(field{1}) = {info.(field{1})};
                end
            else
                % Handle missing fields
                if ~isfield(unique_fields, field{1})
                    unique_fields.(field{1}) = {};
                end
            end
        end

        % Add filename to the list
        if isfield(unique_fields, 'Filename')
            unique_fields.Filename = [unique_fields.Filename, {file_path}];
        else
            unique_fields.Filename = {file_path};
        end
    end

    % Convert all fields' values to unique sets
    fields_names = fieldnames(unique_fields);
    for field = fields_names'
        unique_fields.(field{1}) = unique(unique_fields.(field{1}));
    end

    % Convert the struct to a table
    dicom_fields = struct2table(unique_fields, 'AsArray', true);
end
