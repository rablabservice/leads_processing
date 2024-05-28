function [ref_regions, suvr_files] = get_suvr_files(pet_dir, ref_region_file)
    % Return reference regions and SUVR filenames for a PET scan
    %
    % Parameters
    % ----------
    % pet_dir : char or str
    %     Path to the directory containing the PET scan data
    % ref_region_file : char or str
    %     Path to the CSV file containing reference regions
    %
    % Returns
    % -------
    % ref_regions : table
    %     Table containing reference regions used to make SUVRs for the
    %     PET tracer
    % suvr_files : cell
    %     Cell array of SUVR filenames
    % ------------------------------------------------------------------
    arguments
        pet_dir {mustBeFolder}
        ref_region_file {mustBeText} = ''
    end

    % Format paths
    pet_tag = get_scan_tag(pet_dir);
    [~, tracer] = parse_scan_tag(pet_tag);

    % Get the list of reference regions, and filter it by tracer
    if isempty(ref_region_file)
        code_dir = fileparts(fileparts(mfilename('fullpath')));
        ref_region_file = fullfile(code_dir, 'config', 'ref_regions.csv');
    end
    ref_regions = readtable(ref_region_file);
    ref_regions = ref_regions(strcmp(ref_regions.tracer, tracer), :);

    % Create a struct to hold the SUVR filenames
    suvr_fields = cellfun( ...
        @(x) append('suvr_', strrep(x, '-', '_')), ...
        ref_regions.ref_region, ...
        'UniformOutput', false ...
    );

    suvr_files = cell2struct( ...
        arrayfun( ...
            @(x) fullfile(pet_dir, append('r', pet_tag, '_suvr-', x, '.nii')), ...
            ref_regions.ref_region, ...
            'UniformOutput', false ...
        ), ...
    suvr_fields);

    % Ensure each value of suvr_files is text, not a cell array
    for k = 1:numel(suvr_fields)
        suvr_files.(suvr_fields{k}) = suvr_files.(suvr_fields{k}){1};
    end
end
