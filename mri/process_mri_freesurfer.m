function process_mri_freesurfer(raw_mrif, mri_dir, segment_brainstem, overwrite)
    % Process raw MRIs using FreeSurfer recon-all and segmendBS.sh
    %
    % Parameters
    % ----------
    % raw_mrif : char or str array
    %     Full path to the raw MRI file
    % mri_dir : char or str array
    %     Full path to the FreeSurfer directory that will be created
    % segment_brainstem : logical
    %     If true, segment the brainstem using segmentBS.sh
    % overwrite : logical
    %     If true, overwrite the FreeSurfer directory if it already exists
    % ------------------------------------------------------------------
    arguments
        raw_mrif {mustBeFile}
        mri_dir {mustBeFolder}
        segment_brainstem logical = true
        overwrite logical = false
    end

    % Format paths
    fs_version = 'freesurfer_7p1';
    raw_mrif = abspath(raw_mrif);
    mri_dir = abspath(mri_dir);
    fs_dir = fullfile(mri_dir, fs_version);

    % Check if the FreeSurfer directory already exists
    if exist(fs_dir, 'dir')
        if overwrite
            fprintf('- Removing existing FreeSurfer directory: %s\n', fs_dir);
            rmdir(fs_dir, 's');
        else
            fprintf('- FreeSurfer directory exists, will not rerun\n');
            return
        end
    end

    % Run recon-all
    cmd_fs=char(append('recon-all -all -i ', raw_mrif, ' -sd ', mri_dir, ' -s ', fs_version));
    fprintf('- Processing MRI through FreeSurfer\n')
    fprintf('  $ %s\n', cmd_fs);
    run_system_cmd(cmd_fs);

    % Segment the brainstem
    if segment_brainstem
        % Run segmentBS.sh
        cmd_bs=char(append('segmentBS.sh ', fs_version, ' ', mri_dir));
        fprintf('  $ %s\n', cmd_bs);
        run_system_cmd(cmd_bs);
    end
end
