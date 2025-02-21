function all_exist = processed_wmap_exist(run_dir)
    % Return true if W-map exists in the directory
    %
    % Parameters
    % ----------
    % run_dir : char or str
    %     The directory that contains the processed MRI or PET data
    %
    % Returns
    % -------
    % all_exist : logical
    %    true if all processed MRI files are present, else false
    % ------------------------------------------------------------------
    arguments
        run_dir {mustBeText} = ''
    end
    run_dir = abspath(run_dir);

    suvr_files = {};
    wmap_files = {};
    
    code_dir = fileparts(fileparts(mfilename('fullpath')));
    ref_region_file = fullfile(code_dir, 'config', 'ref_regions.csv');
    ref_regions = readtable(ref_region_file);
    
    % Get scan info
    scan_tag = get_scan_tag(run_dir);
    [~, tracer] = parse_scan_tag(scan_tag);

    % Get W-map files
    if ~strcmp(tracer,'MRI-T1')
        ref_region = ref_regions(strcmp(ref_regions.tracer, tracer), :);
        ref_region = ref_region.ref_region(1); % cross-sectional ref region selection
        wmap_file = fullfile(run_dir, append('W-map_wr', scan_tag, '_suvr-', ref_region, '.nii')); 
    else
        wmap_file = fullfile(run_dir, append('W-map_s8mwc1', scan_tag, '_nu.nii')); 
    end

    all_exist = all(isfile(cellvec(wmap_file)));
