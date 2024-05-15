function setup_leads_processing(data_dir, overwrite, process_unused_mris, cleanup)
    % High-level function that sets up the LEADS processing pipeline
    % ------------------------------------------------------------------
    arguments
        data_dir {mustBeFolder} = '/mnt/coredata/processing/leads/data'
        overwrite logical = false
        process_unused_mris logical = false
        cleanup logical = true
    end

    % Format paths to the data directories
    data_dir = abspath(data_dir);
    newdata_dir = fullfile(data_dir, 'newdata');
    raw_dir = fullfile(data_dir, 'raw');

    % Set path to the python interpreter
    python = "/home/mac/dschonhaut/mambaforge/envs/nipy311/bin/python";
    pyenv(Version=python);

    % Add paths to the Python environment
    code_dir = fileparts(mfilename('fullpath'));
    if count(py.sys.path, code_dir) == 0
        insert(py.sys.path, int32(0), code_dir);
    end

    % Print the module header
    title = 'SETUP MODULE';
    subtitle = 'Prepare new MRI and PET scans for processing';
    print_title(title, subtitle);

    % Unzip files in newdata
    zipfs = glob_sort_mtime(fullfile(newdata_dir, '*.zip'));
    for ii = 1:numel(zipfs)
        fprintf('- Unzipping %s\n', zipfs{ii});
        unzip(zipfs{ii}, newdata_dir);
        % Remove the zip file
        delete(zipfs{ii});
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

    % Exit early
    return

    % % Save CSVs files of MRI and PET scans in the raw directory, and
    % % indicate which scans are scheduled for processing
    % cmd = append(python, ' ', fullfile(code_dir, 'select_scans_to_process.py'));
    % if overwrite
    %     cmd = append(cmd, ' -o');
    % end
    % if process_unused_mris
    %     cmd = append(cmd, ' -p');
    % end
    % run_system_cmd(cmd)

    % Create processed scan directories for MRI and PET scans that need
    % to be processed, link each PET scan to its closest MRI, and copy
    % PET niftis from their raw to processed directories
    cmd = append(python, ' ', fullfile(code_dir, 'setup_processed_scan_dirs.py'));
    if overwrite
        cmd = append(cmd, ' -o');
    end
    run_system_cmd(cmd);

    % Print the module footer
    print_footer('Setup module complete');
end
