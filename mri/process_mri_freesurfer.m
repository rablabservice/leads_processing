function process_mri_freesurfer(raw_mrif, mri_dir, segment_brainstem, overwrite, verbose)
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
    % verbose : logical
    %     If true, print diagnostic information
    % ------------------------------------------------------------------
    arguments
        raw_mrif {mustBeText}
        mri_dir {mustBeText}
        segment_brainstem logical = true
        overwrite logical = false
        verbose logical = true
    end

    % Format inputs
    fs_version = 'freesurfer_7p1';
    raw_mrif = abspath(raw_mrif);
    mri_dir = abspath(mri_dir);
    fs_dir = fullfile(mri_dir, fs_version);

    % Check if the FreeSurfer directory already exists
    if exist(fs_dir, 'dir')
        if overwrite
            rmdir(fs_dir, 's');
        else
            error('%s already exists. Set overwrite to true to delete it.', fs_dir);
        end
    end

    % First run recon-all
    cmd_fs=char(append('recon-all -all -i ', raw_mrif, ' -sd ', mri_dir, ' -s ', fs_version));
    if verbose
        disp(cmd_fs);
    end
    system(cmd_fs);

    % Then run segmentBS.sh
    if segment_brainstem
        cmd_bs=char(append('segmentBS.sh ', fs_version, ' ', mri_dir));
        if verbose
            disp(cmd_bs);
        end
        system(cmd_bs);
    end
end
