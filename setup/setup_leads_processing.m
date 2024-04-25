function setup_leads_processing(data_dir, overwrite, cleanup, verbose)
    % High-level function that sets up the LEADS processing pipeline
    % ------------------------------------------------------------------
    arguments
        data_dir {mustBeFolder} = '/mnt/coredata/processing/leads/data'
        overwrite logical = false
        cleanup logical = true
        verbose logical = true
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

    % Convert dicoms to nifti
    convert_dicoms(newdata_dir, verbose);

    % Move data from newdata to raw
    py.move_newdata_to_raw.move_newdata_to_raw(...
        newdata_dir=newdata_dir, ...
        raw_dir=raw_dir, ...
        overwrite=overwrite, ...
        cleanup=cleanup, ...
        verbose=verbose ...
    );

    % % Find scans that need to be processed
    % py.find_scans_to_process.find_scans_to_process(...
    %     raw_dir=raw_dir, ...
    %     processed_dir=processed_dir, ...
    %     overwrite=overwrite, ...
    %     verbose=verbose ...
    % );

    % Create processed PET directories, link them to their appropriate
    % MRI directories, and copy in the raw PET niftis
    cmd = append(python, ' ', fullfile(code_dir, 'create_pet_proc_dirs.py'));
    if overwrite
        cmd = append(cmd, ' -o');
    end
    if ~verbose
        cmd = append(cmd, ' -q');
    end
    run_system_cmd(cmd);
end
