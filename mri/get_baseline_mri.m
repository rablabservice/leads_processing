function baseline_nuf = get_baseline_mri(subj_dir, scan_type)
    % Given path to subject dir, return full path to the baseline nu.nii
    % ------------------------------------------------------------------
    arguments
        subj_dir {mustBeFolder}
        scan_type {mustBeText} = 'MRI-T1'
    end

    % Format inputs
    subj_dir = char(abspath(subj_dir));
    scan_type = char(scan_type);

    % Find all processed scan directories
    files = dir(subj_dir);
    dirs = {files([files.isdir]).name};
    dirs = dirs(~ismember(dirs, {'.', '..'}));

    % Filter directories starting with the MRI tag
    scan_dirs = dirs(startsWith(dirs, scan_type));

    % Extract MRI dates
    scan_dates = cellfun(@(x) strsplit(x, '_'), scan_dirs, 'UniformOutput', false);
    scan_dates = cellfun(@(x) x{2}, scan_dates, 'UniformOutput', false);

    % Find the baseline MRI
    date_nums = cellfun(@(x) datenum(x, 'yyyy-mm-dd'), scan_dates);
    baseline_scan_date = datestr(min(date_nums), 'yyyy-mm-dd');
    baseline_scan_dir = fullfile(subj_dir, strjoin({scan_type, baseline_scan_date}, '_'));
    scan_tag = get_scan_tag(baseline_scan_dir);
    baseline_nuf = fullfile(baseline_scan_dir, strjoin({scan_tag, 'nu.nii'}, '_'));
    baseline_nuf = abspath(baseline_nuf);
end
