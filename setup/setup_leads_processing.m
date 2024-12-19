function setup_leads_processing( ...
    proj_dir, ...
    skip_newdata, ...
    wipe_newdata, ...
    setup_processed, ...
    dry_run, ...
    schedule_all_mris, ...
    overwrite_raw, ...
    schedule_overwrite, ...
    overwrite_processed ...
)
    % Setup the LEADS processing pipeline.
    %
    % Running this function is the precursor to processing MRIs through
    % FreeSurfer and SPM (process_mris.m) and to processing PET scans
    % through the MRI-based pipeline (process_pets.m). See below for an
    % overview everything that happens during setup.
    %
    % Overview
    % --------
    % Part I (newdata -> raw)
    % 1.  Unzip any zip files in newdata
    % 2.  Convert newdata dicoms to nifti
    % 3.  Move MRI and PET scan directories from newdata to raw
    %
    % Part II (raw -> processed)
    % 4.  Save CSV files listing all MRI and PET scans in raw,
    %     respectively. Indicate which scans need to be processed
    %     (column 'scheduled_for_processing'). Once these CSV files are
    %     saved (in proj_dir/metadata/scans_to_process), the user can
    %     manually edit them for greater control over which scans will
    %     be processed. process_mris.m and process_pets.m will read
    %     these files to determine which scans to process.
    % 5.  Create processed scan directories for MRI and PET scans that
    %     are scheduled for processing.
    %     I.   Each processed scan directory contains a symlink called
    %          'raw' that points back to the relevant raw data directory
    %     II.  Each PET scan is also linked to its closest MRI scan,
    %          to which it will later be coregistered
    % ------------------------------------------------------------------
    arguments
        proj_dir {mustBeFolder} = '/mnt/coredata/processing/leads'
        skip_newdata logical = false
        wipe_newdata logical = true
        setup_processed logical = true
        dry_run logical = false
        schedule_all_mris logical = false
        overwrite_raw logical = false
        schedule_overwrite logical = false
        overwrite_processed logical = false
    end

    % Print the module header
    fprintf('\n--------------------------------------------\n');
    fprintf('START SETUP MODULE\n');
    fprintf('Prepare new MRI and PET scans for processing\n')
    fprintf('--------------------------------------------\n\n');

    % Format paths to project directories
    proj_dir = abspath(proj_dir);
    newdata_dir = fullfile(proj_dir, 'data', 'newdata');
    raw_dir = fullfile(proj_dir, 'data', 'raw');
    scans_to_process_dir = fullfile(proj_dir, 'metadata', 'scans_to_process');
    code_dir = fullfile(proj_dir, 'code');

    % Get path to the Python interpreter
    python = '/mnt/coredata/Projects/Resources/dscode/bin/python';

    % Part I: newdata -> raw
    % Unzip any zip files in newdata, convert dicoms to nifti, and move
    % scan directories from newdata to raw
    if ~skip_newdata
        newdatafs = dir(newdata_dir);
        newdatafs = newdatafs(~startsWith({newdatafs.name}, '.'));
        if isempty(newdatafs)
            fprintf('- %s is empty, skipping ahead\n', newdata_dir);
        else
            % Unzip newdata files
            cmd = append(python, fullfile(code_dir, 'setup', 'unzip_files_in_dir.py'));
            cmd = append(cmd, ' -d ', newdata_dir);
            fprintf('- Unzipping zip files\n');
            fprintf('  (unzip_files_in_dir.py)\n');
            system(cmd);

            % Convert dicoms to nifti
            fprintf('- Converting newdata dicoms to nifti\n');
            fprintf('  (convert_dicoms.m)\n');
            convert_dicoms(newdata_dir);

            % Move scans from newdata to raw
            cmd = append(python, fullfile(code_dir, 'setup', 'move_newdata_to_raw.py'));
            cmd = append(cmd, ' --newdata ', newdata_dir);
            cmd = append(cmd, ' --raw ', raw_dir);
            if overwrite_raw
                cmd = append(cmd, ' -o');
            end
            if ~wipe_newdata
                cmd = append(cmd, ' --no-clean');
            end
            fprintf('- Moving scans from ''newdata'' to ''raw''\n');
            fprintf('  (move_newdata_to_raw.py)\n');
            system(cmd);
        end
    end

    % Part II: raw -> processed
    % Decide which MRI and PET scans we'll process, then create
    % processed scan directories for them

    % Save CSVs files of MRI and PET scans in the raw directory, and
    % indicate which scans are scheduled for processing
    cmd = append(python, fullfile(code_dir, 'setup', 'select_scans_to_process.py'));
    cmd = append(cmd, ' -p ', proj_dir);
    if schedule_all_mris
        cmd = append(cmd, ' -a');
    end
    if schedule_overwrite
        cmd = append(cmd, ' -o');
    end
    if dry_run
        cmd = append(cmd, ' --no-save-csv');
    end
    if dry_run
        fprintf('- Scheduling scans to process (dry run, no CSVs will be saved)\n');
    else
        fprintf('- Scheduling scans to process\n');
    end
    fprintf('  (select_scans_to_process.py)\n');
    system(cmd);

    % Create processed scan directories for MRI and PET scans that need
    % to be processed, link each PET scan to its closest MRI, and copy
    % PET niftis from their raw to processed directories
    if setup_processed && ~dry_run
        cmd = append(python, fullfile(code_dir, 'setup', 'make_processed_scan_dirs.py'));
        cmd = append(cmd, ' --scans_to_process_dir ', scans_to_process_dir);
        if overwrite_processed
            cmd = append(cmd, ' -o');
        end
        fprintf('- Setting up scan directories ahead of processing\n');
        fprintf('  (make_processed_scan_dirs.py)\n');
        system(cmd);
    end

    % Print the module footer
    fprintf('\n----------------\n');
    fprintf('END SETUP MODULE\n\n');
end
