function setup_leads_processing(data_dir, overwrite, cleanup)
    % High-level function that sets up the LEADS processing pipeline
    % ------------------------------------------------------------------
    arguments
        data_dir {mustBeFolder} = '/mnt/coredata/processing/leads/data'
        overwrite logical = false
        cleanup logical = true
    end

    % Format paths to the data directories
    data_dir = abspath(data_dir);
    newdata_dir = fullfile(data_dir, 'newdata');
    raw_dir = fullfile(data_dir, 'raw');
    processed_dir = fullfile(data_dir, 'processed');

    % Set path to the python interpreter
    python = "/home/mac/dschonhaut/mambaforge/envs/nipy311/bin/python";
    pyenv(Version=python);

    % Add paths to the Python environment
    code_dir = fileparts(mfilename('fullpath'));
    if count(py.sys.path, code_dir) == 0
        insert(py.sys.path, int32(0), code_dir);
    end

    % Print welcome message
    title = 'SETUP MODULE';
    subcap = 'Prepare new MRI and PET scans for processing';
    border_len = max(length(title), length(subcap));
    border_char = '~';
    fprintf('%s\n', repmat(border_char, 1, border_len));
    fprintf('%s\n\n', repmat(border_char, 1, border_len));
    fprintf('%s\n\n', title);
    fprintf('%s\n\n', subcap);
    while border_len > 0
        fprintf('%s\n', repmat(border_char, 1, border_len));
        border_len = border_len - 4;
    end

    % Convert dicoms to nifti
    convert_dicoms(newdata_dir);

    % Move data from newdata to raw
    py.move_newdata_to_raw.move_newdata_to_raw(...
        newdata_dir=newdata_dir, ...
        raw_dir=raw_dir, ...
        overwrite=overwrite, ...
        cleanup=cleanup ...
    );

    % % Save a CSV file of MRI + PET scans that need to be processed
    % cmd = append(python, ' ', fullfile(code_dir, 'find_scans_to_process.py'));
    % run_system_cmd(cmd)

    % % Create processed scan directories for MRI and PET scans that need
    % % to be processed, link each PET scan to its closest MRI, and copy
    % % PET niftis from their raw to processed directories
    % cmd = append(python, ' ', fullfile(code_dir, 'create_pet_proc_dirs.py'));
    % if overwrite
    %     cmd = append(cmd, ' -o');
    % end
    % run_system_cmd(cmd);
end
