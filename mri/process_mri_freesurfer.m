function outfiles = process_mri_freesurfer( ...
    raw_mrif, ...
    mri_dir, ...
    fid, ...
    overwrite, ...
    segment_brainstem ...
)
    % Process raw MRIs using FreeSurfer recon-all and segmentBS.sh
    %
    % Overview
    % --------
    % 1.  Run recon-all on the raw MRI file
    % 2.  Optionally segment the brainstem using segmentBS.sh
    %
    % Parameters
    % ----------
    % raw_mrif : char or str array
    %     Full path to the raw MRI file
    % mri_dir : char or str array
    %     Full path to the FreeSurfer directory that will be created
    % fid : int, optional
    %     File identifier for logging (default is 1 for stdout)
    % overwrite : logical
    %     If true, overwrite the FreeSurfer directory if it already exists
    % segment_brainstem : logical
    %     If true, segment the brainstem using segmentBS.sh
    % ------------------------------------------------------------------
    arguments
        raw_mrif {mustBeFile}
        mri_dir {mustBeFolder}
        fid {mustBeNumeric} = 1
        overwrite logical = false
        segment_brainstem logical = true
    end

    % Format paths
    fs_version = 'freesurfer_7p1';
    raw_mrif = abspath(raw_mrif);
    mri_dir = abspath(mri_dir);
    fs_dir = fullfile(mri_dir, fs_version);
    fs_link = fullfile(mri_dir, 'freesurfer');

    % If processing is already complete and overwrite is false, get
    % the struct of FreeSurfer *.mgz files that we care about and return
    if freesurfer_files_exist(mri_dir)
        if overwrite
            log_append(fid, sprintf('- Removing existing FreeSurfer directory: %s', fs_dir));
            rmdir(fs_dir, 's');
        else
            log_append(fid, '- FreeSurfer processing already complete, will not rerun');
            outfiles = get_freesurfer_files(mri_dir, 'mgz');
            return
        end
    end

    % Run recon-all
    cmd_fs=char(append('recon-all -all -i ', raw_mrif, ' -sd ', mri_dir, ' -s ', fs_version));
    log_append(fid, '- Processing MRI through FreeSurfer');
    log_append(fid, sprintf('    $ %s', cmd_fs));
    system(cmd_fs);

    % Symlink to the freesurfer directory
    if isfolder(fs_dir)
        if system(sprintf('test -L %s', fs_link)) == 0
            delete(fs_link);
        end
        system(sprintf('ln -s %s %s', fs_dir, fs_link));
    end

    % Delete the link to fsaverage
    fsaverage_link = fullfile(mri_dir, 'fsaverage');
    if system(sprintf('test -L %s', fsaverage_link)) == 0
        delete(fsaverage_link);
    end

    % Segment the brainstem
    if segment_brainstem
        % Run segmentBS.sh
        cmd_bs=char(append('segmentBS.sh ', fs_version, ' ', mri_dir));
        log_append(fid, '- Segmenting brainstem into subregions');
        log_append(fid, sprintf('    $ %s', cmd_bs));
        system(cmd_bs);
    end
end
