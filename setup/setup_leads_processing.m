function setup_leads_processing( ...
    proj_dir, ...
    overwrite, ...
    cleanup_newdata, ...
    process_all_mris ...
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
    %     III. PET niftis are copied from raw to processed, where they
    %          are renamed with standard naming conventions
    %     IV.  MRI files are not copied from raw to processed, because
    %          FreeSurfer will do this as part of recon-all
    % ------------------------------------------------------------------
    arguments
        proj_dir {mustBeFolder} = '/mnt/coredata/processing/leads'
        overwrite logical = false
        cleanup_newdata logical = true
        process_all_mris logical = false
    end

    % Print the module header
    fprintf('\n--------------------------------------------\n');
    fprintf('START SETUP MODULE\n');
    fprintf('Prepare new MRI and PET scans for processing\n')
    fprintf('--------------------------------------------\n\n');

    % Format paths to project directories
    proj_dir = abspath(proj_dir);
    data_dir = fullfile(proj_dir, 'data');
    newdata_dir = fullfile(data_dir, 'newdata');
    raw_dir = fullfile(data_dir, 'raw');
    scans_to_process_dir = fullfile(proj_dir, 'metadata', 'scans_to_process');

    % Get path to the Python interpreter
    python = '/home/mac/dschonhaut/mambaforge/envs/nipy311/bin/python ';

    % Note that Python scripts that we need to run must be in the same
    % directory as this calling script
    code_dir = fileparts(mfilename('fullpath'));

    % Part I: newdata -> raw
    % Unzip any zip files in newdata, convert dicoms to nifti, and move
    % scan directories from newdata to raw
    newdatafs = dir(newdata_dir);
    newdatafs = newdatafs(~startsWith({newdatafs.name}, '.'));
    if isempty(newdatafs)
        fprintf('- %s is empty, skipping ahead\n', newdata_dir);
    else
        % Unzip newdata files
        cmd = append(python, fullfile(code_dir, 'unzip_files_in_dir.py'));
        cmd = append(cmd, ' -d ', newdata_dir);
        system(cmd);

        % Convert dicoms to nifti
        fprintf('- Converting newdata dicoms to nifti\n');
        convert_dicoms(newdata_dir);

        % Move scans from newdata to raw
        cmd = append(python, fullfile(code_dir, 'move_newdata_to_raw.py'));
        cmd = append(cmd, ' --newdata ', newdata_dir);
        cmd = append(cmd, ' --raw ', raw_dir);
        if overwrite
            cmd = append(cmd, ' -o');
        end
        if ~cleanup_newdata
            cmd = append(cmd, ' --no-clean');
        end
        fprintf('  $ %s', cmd);
        system(cmd);
    end

    % Part II: raw -> processed
    % Decide which MRI and PET scans we'll process, then create
    % processed scan directories for them

    % Save CSVs files of MRI and PET scans in the raw directory, and
    % indicate which scans are scheduled for processing
    cmd = append(python, fullfile(code_dir, 'select_scans_to_process.py'));
    if process_all_mris
        cmd = append(cmd, ' -a');
    end
    if overwrite
        cmd = append(cmd, ' -o');
    end
    fprintf('- Selecting scans to process\n');
    system(cmd);

    % Create processed scan directories for MRI and PET scans that need
    % to be processed, link each PET scan to its closest MRI, and copy
    % PET niftis from their raw to processed directories
    cmd = append(python, fullfile(code_dir, 'make_processed_scan_dirs.py'));
    cmd = append(cmd, ' --scans_to_process_dir ', scans_to_process_dir);
    if overwrite
        cmd = append(cmd, ' -o');
    end
    fprintf('- Setting up scan directories ahead of processing\n');
    system(cmd);

    % Print the module footer
    fprintf('\n----------------\n');
    fprintf('END SETUP MODULE\n\n');

end
