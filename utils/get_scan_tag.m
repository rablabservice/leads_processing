function scan_tag = get_scan_tag(infile)
    % Given path to a scan processing directory or file directly within
    % it, return the '<subj>_<scan_type>_<scan_date>' scan tag
    % ------------------------------------------------------------------
    arguments
        infile {mustBeText}
    end
    infile = abspath(infile);

    % Get the processed scan directory from infile
    [~, name, ext] = fileparts(infile);
    if contains([name ext], '.')
        scan_proc_dir = fileparts(infile);
    else
        scan_proc_dir = fileparts([infile '/']);
    end

    % Split the path into parts
    parts = strsplit(scan_proc_dir, filesep);
    parts = parts(~cellfun('isempty', parts));
    subj = char(parts{end-1});
    parts = parts{end};
    parts = strsplit(parts, '_');
    scan_type = char(parts{1});
    scan_date = char(parts{2});
    scan_tag = strjoin({subj, scan_type, scan_date}, '_');
end
